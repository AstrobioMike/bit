"""
The check made up front by anything that is about to download something: is there an internet connection at all,
and do the sites this run downloads from answer?

require_internet() is given the URLs about to be downloaded from (the same constants the downloads use, so the two
can't drift apart) and stops with a message when they can't be reached:

  - nothing answers at all (neither those sites nor a couple of large ones that are always up): no connection;
  - the internet is reachable but some of those sites aren't, also when given longer to answer: they are named.

When reference databases are missing, both messages also list them and how to get each on another machine.

A site counts as reachable when it gives any HTTP response at all, including an error status (403, 404, 405, 500):
only no answer (a timeout, the name not resolving, the connection being refused or reset) means it can't be reached.
The requests go through urllib, so the proxy environment variables are respected. Nothing is checked unless a
download is about to start, so runs that have everything they need work without a connection.
"""

import http.client
import ssl
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from bit.modules.general import notify_premature_exit, report_message

# large sites that are always up, asked alongside the ones a run needs: if none of them answers either, there is no
# connection (rather than a site being down)
GENERAL_SITES = ("https://github.com", "https://www.google.com")

TIMEOUT = 5         # seconds a site is given to answer
SLOW_TIMEOUT = 30   # what a site that didn't is given on a second try, before it is said not to be answering

NO_CONNECTION = ("This would require downloading something from the internet, but it seems there is no connection "
                 "right now. So we can't proceed :(")


def site(url):
    """ the site a URL is on, as scheme://host (a bare host name is taken as https) """

    parts = urllib.parse.urlsplit(url if "://" in url else f"https://{url}")
    return f"{parts.scheme}://{parts.netloc}"


def host(url):
    """ the host name of a URL, as it is named in the messages """

    return urllib.parse.urlsplit(site(url)).hostname


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """ a redirect is an answer from the site asked, so it isn't followed (urllib then raises it as an HTTPError) """

    def redirect_request(self, *args, **kwargs):
        return None


def probe(url, timeout=TIMEOUT):
    """
    Whether the site answers: True for any HTTP response, whatever its status. HEAD is tried first, then GET if the
    connection was made but the HEAD request itself got no proper response.
    """

    opener = urllib.request.build_opener(_NoRedirect)
    for method in ("HEAD", "GET"):
        request = urllib.request.Request(site(url) + "/", method=method, headers={"User-Agent": "bit"})
        try:
            with opener.open(request, timeout=timeout):
                return True
        except urllib.error.HTTPError:
            return True
        except urllib.error.URLError as problem:
            # the connection couldn't be made (name not resolved, refused, timed out); a certificate that isn't
            # accepted is still an answer from the site, and the download itself will say what is wrong with it
            return isinstance(problem.reason, ssl.SSLCertVerificationError)
        except (http.client.HTTPException, ConnectionError):
            continue    # connected, but the request was dropped: try the next method
        except (OSError, ValueError):
            return False
    return False


def check_sites(urls, timeout=TIMEOUT):
    """ {site: whether it answers} for the sites of the given URLs, asked all at once """

    sites = list(dict.fromkeys(site(url) for url in urls))
    if not sites:
        return {}
    with ThreadPoolExecutor(max_workers=len(sites)) as pool:
        return dict(zip(sites, pool.map(lambda s: probe(s, timeout), sites)))


def _listed(names):
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def require_internet(urls, what="the data", missing_dbs=None, ref_dbs_dir=None, ref_dbs_in_output_dir=False,
                     for_a_run=True, env_note=None):
    """
    Stops with a message unless the sites of the given URLs (what is about to be downloaded from) can be reached.

    what:         what is being downloaded, for the closing line ("the db", "the data"); None when it isn't something
                  that could be fetched on another machine and moved here (e.g. conda environments)
    missing_dbs:  names of the reference databases (bit.modules.ref_dbs) that aren't set up yet under ref_dbs_dir;
                  the message then lists them and how to get each on another machine
    env_note:     lines saying how the conda environments still to be created could be set up another way (when
                  some are), printed at the end of either message
    """

    needed = list(dict.fromkeys(site(url) for url in urls))
    if not needed:
        return

    # all asked at once, so this takes as long as the slowest to answer
    answers = check_sites(needed + list(GENERAL_SITES))
    down = [s for s in needed if not answers[s]]
    if not down:
        return

    # with a connection, a site that hasn't answered may just be slow (the downloads themselves wait far longer
    # than this check does): it is asked once more, with more time, before the run is stopped over it
    if any(answers.values()):
        answers.update(check_sites(down, SLOW_TIMEOUT))
        down = [s for s in needed if not answers[s]]
        if not down:
            return

    db_note = []
    if missing_dbs:
        from bit.modules.ref_dbs import missing_dbs_note
        db_note = missing_dbs_note(ref_dbs_dir, missing_dbs, in_output_dir=ref_dbs_in_output_dir, for_a_run=for_a_run)

    def say(text, color="yellow"):
        report_message(text, color=color, initial_indent="    ", subsequent_indent="    ")

    if not any(answers.values()):
        say(NO_CONNECTION)
        if db_note:
            print("\n" + "\n".join(db_note))
        if env_note:
            print("\n" + "\n".join(env_note))
    else:
        names = [host(s) for s in down]
        one = len(names) == 1
        say(f"The internet is reachable, but {_listed(names)} {'isn' if one else 'aren'}'t answering right now; "
            f"{'it' if one else 'they'} may be down, or blocked from this network.")
        if db_note:
            print("\n" + "\n".join(db_note))
        if what:
            say(f"You can try again later, or you can get {what} from another machine that can access the site(s) "
                "and then move the data here.", color=None)
        else:
            say("You can try again later, or run this from a machine or network that can access the site(s).",
                color=None)
        if env_note:
            print("\n" + "\n".join(env_note))

    notify_premature_exit()


def note_for_a_dry_run(urls):
    """ for a dry run, which checks nothing: says which sites the real run would download from """

    names = list(dict.fromkeys(host(url) for url in urls))
    if names:
        print(f"    The run would need the internet, to download from: {_listed(names)}")


def known_path_env_note(env_name, yaml_path, env_path):
    """ for a conda environment bit creates itself at a known path (e.g. bit assemble's): where it would need to go """

    return [
        f"    The conda environment this needs ({env_name}) isn't set up yet. If you can get it there yourself, you could",
        f"    use or match its yaml ({yaml_path}) and install it here: {env_path}",
    ]


def shared_conda_prefix_note(env_names, create_only_hint=None):
    """
    For the conda environments snakemake creates (it names their folders itself, so one installed by hand wouldn't be
    found): creating them from a machine that shares this filesystem and can reach the channels.
    """

    how = create_only_hint or "running the same command there creates them"
    return [
        f"    The conda environments this needs ({', '.join(env_names)}) aren't set up yet. If this machine shares its",
        f"    filesystem with one that can access the site(s), {how}.",
    ]
