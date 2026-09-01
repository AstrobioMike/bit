import os
import sys
import time
import gzip
import random
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests # type: ignore
from dataclasses import dataclass
from tqdm import tqdm # type: ignore
from bit.modules.general import (color_text, check_files_are_found, wprint,
                                 attempt_to_make_dir, report_message)
from bit.modules.taxonomy.tax_select import (TaxonNotFound, AmbiguousTaxon,
                                             CrossDomainTaxon)
from bit.modules.taxonomy.target_taxon import (resolve_target_taxon_accessions,
                                               expand_target_taxa,
                                               describe_all_expansion,
                                               ensure_source_data,
                                               TargetTaxonError)
from bit.modules.ncbi.parse_ncbi_assembly_summary import parse_ncbi_assembly_summary
from bit.modules.ncbi.get_ncbi_assembly_data import get_ncbi_assembly_data, ncbi_data_table_path


TRANSIENT_STATUS = {429, 500, 502, 503, 504}

max_threads=20
max_retries=10

# statuses that mean we are specifically being rate-limited, as opposed to the
# broader TRANSIENT_STATUS (which also covers plain server hiccups). These get true
# exponential backoff; everything else transient gets a limited backoff
THROTTLE_STATUS = {429}

# ceiling on any single retry sleep (seconds), applied to the throttled path.
MAX_BACKOFF = 30

# separate, larger ceiling for a server-specified Retry-After
MAX_RETRY_AFTER = 300

SAWTOOTH_CYCLE = 5

def dl_ncbi_assemblies(args):

    preflight_checks(args)

    wanted_accs, selections, expansion_note = resolve_targets(args)

    if getattr(args, "dry_run", False):
        report_selection(wanted_accs, selections, expansion_note, args)
        return

    run_data = setup(args, wanted_accs)

    if selections:
        report_selection(wanted_accs, selections, expansion_note, args)

    run_data = parse_main_assembly_table(run_data)

    summarize_search(run_data)

    run_data = download_assemblies(run_data)

    report_finish(run_data)


def preflight_checks(args):

    accessions_file = getattr(args, "wanted_accessions", None)
    target_taxa = getattr(args, "target_taxon", None)

    if not accessions_file and not target_taxa:
        print("")
        wprint(color_text("Nothing to download. Provide a target taxon with `-t`, an accessions file with `-w`, "
                          "or both.", "yellow"))
        print("")
        sys.exit(1)

    if accessions_file:
        check_files_are_found([accessions_file])

    # a dry run reports and stops, so it must not create anything
    if not getattr(args, "dry_run", False):
        if args.output_dir and not os.path.exists(args.output_dir):
            attempt_to_make_dir(args.output_dir)

    get_ncbi_assembly_data(quiet=True)

    if target_taxa:
        # fetches the selection asset (and, for GTDB, the NCBI table it screens against)
        ensure_source_data(getattr(args, "source", "ncbi"))


def _selection_kwargs(args):
    """The selection knobs shared by every `-t` target in a run."""
    source = getattr(args, "source", "ncbi")

    # Deliberately NOT the source's own default (`None`), which would make a GTDB pull
    # species-representatives-only. Two reasons:
    #   1. --derep-rank is `auto` here, so volume is already controlled by dereplication.
    #      Reps-only was a volume brake, and on top of derep it only narrows what each
    #      group can be represented BY -- quietly excluding a possibly-better genome for
    #      no benefit.
    #   2. `bit get-accs-from-gtdb` already defaults to all genomes, so deferring to the
    #      source default here would make this subcommand the odd one out.
    # With this, --representatives-only means the same thing for both sources: restrict
    # to that source's representative/reference genomes.
    reps_only = bool(getattr(args, "representatives_only", False))

    return {
        "target_rank": getattr(args, "target_rank", None),
        "derep_rank": getattr(args, "derep_rank", "auto"),
        "target_domain": getattr(args, "target_domain", None),
        "ncbi_section": getattr(args, "ncbi_section", "refseq"),
        "assembly_levels": getattr(args, "assembly_level", None),
        "reps_only": reps_only,
        "min_completeness": getattr(args, "min_completeness", None),
        "max_contamination": getattr(args, "max_contamination", None),
        "include_rows": False,
    }


def resolve_targets(args):
    """
    Turn `-w` and/or `-t` into ONE deduplicated accession list.

    Returns (accessions, selections, expansion_note). `selections` is one record per
    resolved `-t` target carrying what it selected and how many of those were NEW
    after deduplication -- with a repeatable `-t`, overlapping taxa mean the per-taxon
    counts will not sum to the total, and that difference is the thing worth showing.
    """
    accessions = []
    seen = set()

    accessions_file = getattr(args, "wanted_accessions", None)
    if accessions_file:
        with open(accessions_file, "r") as f:
            for line in f:
                acc = line.strip()
                if acc and acc not in seen:
                    seen.add(acc)
                    accessions.append(acc)

    num_from_file = len(accessions)

    target_taxa = getattr(args, "target_taxon", None)
    if not target_taxa:
        return accessions, [], None

    source = getattr(args, "source", "ncbi")

    try:
        expanded, domains = expand_target_taxa(source, target_taxa)
    except TargetTaxonError as e:
        _exit_with(str(e))
    expansion_note = describe_all_expansion(source, domains)

    kwargs = _selection_kwargs(args)
    selections = []

    for taxon in expanded:
        try:
            taxon_accs, selection = resolve_target_taxon_accessions(
                source, taxon, **kwargs)
        except AmbiguousTaxon as e:
            _exit_with(f"'{e.taxon}' occurs at more than one rank "
                       f"({', '.join(e.ranks_found)}). Specify which is wanted with "
                       f"`-r`.")
        except CrossDomainTaxon as e:
            _exit_with(f"'{e.taxon}' occurs in more than one domain "
                       f"({', '.join(e.domains_found)}), so pulling on the name alone "
                       f"would mix genomes from different domains. Specify which is "
                       f"wanted with `--target-domain` (e.g. `--target-domain "
                       f"{e.domains_found[0]}`).")
        except TaxonNotFound:
            _exit_with(f"Input taxon '{taxon}' doesn't seem to exist at any rank in "
                       f"the {source.upper()} table :(")
        except (TargetTaxonError, ValueError) as e:
            _exit_with(str(e))

        num_new = 0
        for acc in taxon_accs:
            if acc not in seen:
                seen.add(acc)
                accessions.append(acc)
                num_new += 1

        selections.append(TaxonSelection(
            taxon=taxon,
            canonical=selection.canonical,
            resolved_rank=selection.resolved_rank,
            effective_derep_rank=selection.effective_derep_rank,
            num_selected=len(taxon_accs),
            num_new=num_new,
            warnings=list(selection.warnings)))

    return accessions, selections, (expansion_note, num_from_file)


def _exit_with(message):
    print("")
    wprint(color_text(message, "yellow"))
    print("")
    sys.exit(1)


@dataclass
class TaxonSelection:
    """What one `-t` target resolved to, for the reporting layer."""
    taxon: str = None
    canonical: str = None
    resolved_rank: str = None
    effective_derep_rank: str = None
    num_selected: int = 0
    num_new: int = 0
    warnings: list = None


def report_selection(accessions, selections, expansion_note, args):
    """
    What each `-t` resolved to, and the combined total.

    `num_selected` is what the taxonomy core returned for that taxon; `num_new` is what
    survived deduplication against `-w` and against earlier `-t` targets. They differ
    whenever two taxa overlap, which is exactly the case where a bare per-taxon count
    would mislead. Shared with `--dry-run`, so the numbers can't drift apart.
    """
    if getattr(args, "quiet", False):
        return

    note, num_from_file = (expansion_note if expansion_note else (None, 0))

    print("")
    if note:
        wprint(color_text(note, "yellow"))
        print("")

    if num_from_file:
        wprint(f"{num_from_file:,} accession(s) read from "
               f"{color_text(args.wanted_accessions)}")
        print("")

    for sel in selections:
        derep = (f"dereplicated to one genome per {sel.effective_derep_rank}"
                 if sel.effective_derep_rank else "dereplication off")
        wprint(f"-t {color_text(sel.canonical)} "
               f"(resolved to {sel.resolved_rank}; {derep})")

        line = f"{sel.num_selected:,} genome(s) selected"
        overlap = sel.num_selected - sel.num_new
        if overlap:
            line += f", {sel.num_new:,} new ({overlap:,} already counted)"
        wprint("    " + line)

        for warning in (sel.warnings or []):
            report_message(warning, "yellow")
        print("")

    total = len(accessions)
    if getattr(args, "dry_run", False):
        wprint(color_text(f"Total that would be downloaded: {total:,} genome(s) "
                          f"in {args.format} format.", "green"))
    else:
        wprint(color_text(f"Total to download: {total:,} genome(s)", "green"))
    print("")


def setup(args, wanted_accs=None):

    if wanted_accs is None:
        with open(args.wanted_accessions, "r") as f:
            wanted_accs = [line.strip() for line in f if line.strip()]

    run_data = RunData(wanted_format=args.format,
        num_jobs=args.jobs,
        output_dir=args.output_dir,
        wanted_accs=wanted_accs,
        num_wanted=len(wanted_accs),
        ncbi_sub_table_path=Path(args.output_dir) / "wanted-ncbi-accessions-info.tsv",
        not_found_path=Path(args.output_dir) / "ncbi-accessions-not-found.txt",
        not_downloaded_path=Path(args.output_dir) / "ncbi-accessions-not-downloaded.tsv",
        quiet=getattr(args, "quiet", False),
        from_taxon=bool(getattr(args, "target_taxon", None))
        )

    return run_data


def parse_main_assembly_table(run_data):

    if not run_data.quiet:
        print(color_text(f"\n    Targeting {run_data.num_wanted} accession(s) in {run_data.wanted_format} format...\n", "yellow"))

    assembly_summary_file = ncbi_data_table_path()

    run_data = parse_ncbi_assembly_summary(assembly_summary_file, run_data)

    return run_data


@dataclass
class RunData:
    wanted_format: str = None
    num_jobs: int = 10
    output_dir: str = None
    wanted_accs: list = None
    num_wanted: int = 0
    num_found: int = 0
    num_not_found: int = 0
    num_downloaded: int = 0
    num_skipped: int = 0
    num_not_downloaded: int = 0
    ncbi_sub_table_path: str = None
    not_found_path: str = None
    not_downloaded_path: str = None
    quiet: bool = False
    from_taxon: bool = False

    @property
    def not_found_reason(self):
        """
        Why accessions might be missing depends on where they came from.

        Telling someone who passed `-t` that their input "may be invalid" is wrong and
        confusing: they never typed an accession. On that path we just state the count
        and point at the file, the way GToTree does, rather than speculating about the
        cause or suggesting they refresh anything.
        """
        if self.from_taxon:
            return ""
        return "They may be invalid or suppressed."

    @property
    def none_found_hint(self):
        if self.from_taxon:
            return "None of the target taxa could be found."
        return "This is kinda weird. Are the inputs assembly accessions?"


def summarize_search(summary):

    if summary.num_found != summary.num_wanted:
        if summary.num_found > 0:
            print(color_text(f"{' ' * 34}NOTICE", "orange"))
            reason = f" {summary.not_found_reason}" if summary.not_found_reason else ""
            print(f"      {summary.num_not_found} accession(s) not found at NCBI.{reason}")
            print(f"          See '{summary.not_found_path}'.\n")
        else:
            print(color_text(f"{' ' * 34}NOTICE", "orange"))
            print(f"      None of the {summary.num_wanted} target accession(s) were found at NCBI...")
            print(f"      {summary.none_found_hint}\n")
            os.remove(summary.not_found_path)
            sys.exit(1)

        if not summary.quiet:
            print(f"    Remaining total targets: {summary.num_found}\n")


def sleep_backoff(attempt, resp=None, throttled=False):
    """
    Sleep before a retry. Two policies, chosen by WHY the attempt failed:

    throttled=True  -- the server is explicitly rate-limiting us (HTTP 429, or any
                       response carrying Retry-After). Honor Retry-After if given
                       (capped at MAX_RETRY_AFTER), otherwise true exponential
                       backoff capped at MAX_BACKOFF. Backing off progressively is
                       the whole point here: retrying faster accelerates into the
                       thing that's throttling us and tends to extend the penalty
                       window. Retry-After gets the more generous ceiling because
                       it's a specific instruction rather than a guess -- clamping
                       it to MAX_BACKOFF would just burn retries on requests the
                       server has already said it won't serve yet.

    throttled=False -- a generic transient failure (5xx, connection reset, timeout,
                       truncated/empty body, NCBI error page). These are blips or
                       genuinely dead URLs, and there's no politeness argument for
                       stretching the interval forever, so the wait SAWTOOTHS:
                       1, 2, 4, 8, 16, 1, 2, 4, 8, 16, ... That keeps a straggler
                       from parking a pool thread for hours while still spacing
                       requests out, and we find out a file is dead far sooner.
    """
    if throttled:
        if resp is not None:
            retry_after = resp.headers.get("Retry-After")
            if retry_after:
                try:
                    time.sleep(min(float(retry_after), MAX_RETRY_AFTER))
                    return
                except ValueError:
                    pass
        time.sleep(min((2 ** (attempt - 1)) + random.uniform(0, 1), MAX_BACKOFF))
        return

    # sawtooth: exponential within a short cycle, then start over
    step = (attempt - 1) % SAWTOOTH_CYCLE
    time.sleep((2 ** step) + random.uniform(0, 1))


def download_one(target_link, local_dest, retries=max_retries):

    local_path = Path(local_dest)
    tmp_path = local_path.with_name(local_path.name + ".tmp")

    # atomic writes (below) guarantee any 'final file' completed this run
    if local_path.exists() and local_path.stat().st_size > 0:
        return (local_dest, None, "skipped")

    local_path.parent.mkdir(parents=True, exist_ok=True)

    # clear any stale tmp left behind by a prior interrupted run
    tmp_path.unlink(missing_ok=True)

    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(target_link, stream=True, timeout=60)

            if resp.status_code == 404:
                return (local_dest, "Not available in requested format (404)", "failed")

            if resp.status_code in TRANSIENT_STATUS:
                if attempt == retries:
                    return (local_dest, f"HTTP {resp.status_code} after {retries} attempts", "failed_transient")
                # a 429, or any response that bothered to tell us when to come back,
                # is a throttle -- back off properly. A bare 5xx is just a hiccup.
                throttled = (resp.status_code in THROTTLE_STATUS
                             or resp.headers.get("Retry-After") is not None)
                sleep_backoff(attempt, resp, throttled=throttled)
                continue

            resp.raise_for_status()

            content_type = resp.headers.get("Content-Type", "")
            if "xml" in content_type.lower() or "html" in content_type.lower():
                if attempt == retries:
                    return (local_dest, f"NCBI returned an error page after {retries} attempts", "failed_transient")
                sleep_backoff(attempt)
                continue

            with open(tmp_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1024 * 64):
                    fh.write(chunk)

            if tmp_path.stat().st_size == 0:
                tmp_path.unlink(missing_ok=True)
                if attempt == retries:
                    return (local_dest, "Downloaded file was empty", "failed_transient")
                sleep_backoff(attempt)
                continue

            os.replace(tmp_path, local_path)
            return (local_dest, None, "downloaded")

        except (requests.RequestException, OSError) as e:
            tmp_path.unlink(missing_ok=True)
            if attempt == retries:
                return (local_dest, str(e), "failed_transient")
            sleep_backoff(attempt)


def run_download_pass(targets, run_data, desc="Progress"):
    """
    runs one pooled download pass over a list of (target_link, local_dest) tuples.
    returns (permanent_failures, transient_failures, num_skipped) where each
    failure list holds (target_link, local_dest, error) so transient ones can be retried
    """
    permanent = []
    transient = []
    num_skipped = 0

    link_by_dest = {dest: link for link, dest in targets}

    pool = ThreadPoolExecutor(max_workers=min(run_data.num_jobs, max_threads))
    try:
        futures = {
            pool.submit(download_one, link, dest): dest
            for link, dest in targets
        }

        if run_data.quiet:
            desc_buffer = "    "
            ncols = 78
        else:
            desc_buffer = "      "
            ncols = 70
        with tqdm(total=len(targets), desc=f"{desc_buffer}{desc}", unit=" file", ncols=ncols, smoothing=0.05) as pbar:
            for future in as_completed(futures):
                dest, error, status = future.result()
                if status == "failed":
                    permanent.append((link_by_dest[dest], dest, error))
                elif status == "failed_transient":
                    transient.append((link_by_dest[dest], dest, error))
                elif status == "skipped":
                    num_skipped += 1
                pbar.update(1)
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    finally:
        pool.shutdown(wait=True)

    return permanent, transient, num_skipped


def download_assemblies(run_data):

    if not run_data.quiet:
        print(color_text("    Downloading assemblies...\n", "yellow"))

    targets = []
    with open(run_data.ncbi_sub_table_path, "r") as f:
        header = f.readline().strip().split("\t")
        link_idx = header.index("target_link")
        dest_idx = header.index("local_destination")
        for line in f:
            fields = line.strip().split("\t")
            targets.append((fields[link_idx], fields[dest_idx]))

    permanent, transient, num_skipped = run_download_pass(targets, run_data)

    # second pass on transient-only failures
    if transient:
        retry_targets = [(link, dest) for link, dest, _ in transient]
        print(color_text(f"\n    {len(transient)} file(s) failed with transient error messages, doing another pass", "yellow"))
        print(color_text(f"    to see if we can grab them...\n", "yellow"))

        time.sleep(3)
        retry_permanent, retry_transient, retry_skipped = run_download_pass(
            retry_targets, run_data, desc="Progress"
        )
        num_skipped += retry_skipped
        # anything still failing after the retry is final, regardless of category
        permanent.extend(retry_permanent)
        permanent.extend(retry_transient)

    failed = [(dest, error) for _, dest, error in permanent]

    if failed:
        with open(run_data.not_downloaded_path, "w") as fh:
            fh.write("accession\terror\n")
            for dest, error in failed:
                acc = Path(dest).stem.split(".")[0]
                fh.write(f"{acc}\t{error}\n")
        run_data.num_not_downloaded = len(failed)
    else:
        Path(run_data.not_downloaded_path).unlink(missing_ok=True)

    run_data.num_skipped = num_skipped
    run_data.num_downloaded = len(targets) - len(failed)

    return run_data


def report_finish(run_data):

    skipped_note = ""
    if run_data.num_skipped > 0:
        skipped_note = f" ({run_data.num_skipped} already present, skipped)"

    if run_data.num_downloaded == run_data.num_wanted:
        if not run_data.quiet:
            print(color_text(f"\n    All {run_data.num_wanted} file(s) downloaded successfully!{skipped_note}\n", "green"))

    elif run_data.num_downloaded == run_data.num_found:
        if not run_data.quiet:
            print(color_text(f"\n    All {run_data.num_found} found file(s) downloaded successfully!{skipped_note}\n", "yellow"))

    elif run_data.num_not_downloaded > 0:
        print(color_text(f"\n{' ' * 34}NOTICE", "orange"))
        print(f"      {run_data.num_not_downloaded} file(s) failed to download from NCBI. They may not be available")
        print(f"      in the requested format, or it may have been a transient problem.")
        print(f"          See '{run_data.not_downloaded_path}'.")

        if run_data.num_downloaded > 0:
            if not run_data.quiet:
                print(color_text(f"\n\n    The remaining {run_data.num_downloaded} found file(s) downloaded successfully.{skipped_note}\n", "yellow"))
        else:
            print(color_text(f"\n\n    No files were successfully downloaded...{skipped_note}\n", "orange"))
            sys.exit(1)
