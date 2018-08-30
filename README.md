## Bioinformatics Tools (bit)
Here are some useful bioinformatics one-liners and scripts I use frequently enough that it's been worth it to me to have them instantly available everywhere. 

If you clone or download this repository and add it to your PATH, they will available everywhere for you too. If you're not yet sure what your PATH is and what modifying it is all about, [see here](https://astrobiomike.github.io/bash/modifying_your_path#adding-a-directory-to-your-path) 🙂

[astrobiomike.github.io](https://astrobiomike.github.io/)  
[@AstrobioMike](https://twitter.com/AstrobioMike)

## Commands and usage
Each has a help menu accessible by either entering the command alone or by providing `-h` as the only argument.  

|command|function|
|:-----:|-----|
|[bit-remove-wraps](#remove-line-wraps-from-fasta-file)|Remove line wraps from fasta file|
|[bit-count-bases]()|Count number of total bases in fasta file|
|[bit-count-bases-per-seq]()|Count number of bases per sequence|
|[bit-calc-gc-per-sequence]()|Calculate GC of sequences|
|[bit-calc-gc-sliding-window]()|Calculate rolling GC|
|[bit-extract-seqs-by-coords]()|Extract sequences by coordinates|
|[bit-parse-fasta-by-headers]()|Parse fasta by headers|
|[bit-simplify-fasta-headers]()|Simplify sequence headers|
|[bit-reorder-fasta]()|Reorder fasta by headers|
|[bit-calc]()|Save one step calling the command-line calculator (meh)|


Count number of total bases in fasta file
* Calculator
* bit-calculate-

### Remove line wraps from fasta file
```
$ bit-remove-wraps -h
This script removes line wraps from a fasta file.

Usage:
	 bit-remove-wraps input.fasta > new.fasta

```

### Count number of total bases in fasta file

```
$ bit-count-bases -h
This script returns the total number of bases (or amino acids) in a fasta file.

Usage:
	 count-bases input.fasta
```

### Count number of bases per sequence

```
$ bit-count-bases-per-seq -h
usage: bit-count-bases-per-seq [-h] [-i INPUT_FASTA] [-o OUTPUT_FILE]

This script takes a multifasta as input and returns a tab-delimited file with
two columns, header and number of bases or amino acids, for each sequence.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --output_txt_file OUTPUT_FILE
                        Name of output txt file (default: "Num_bps.txt")

required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Original fasta file
```

### Calculate GC of sequences

```
$ bit-calc-gc-per-sequence -h
usage: bit-calc-gc-per-sequence [-h] [-i INPUT_FASTA] [-o OUTPUT_FILE]

This script takes a nucleotide multifasta and returns a tab-delimited file
with 3 columns: header, sequence length, and GC.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --output_txt_file OUTPUT_FILE
                        Name of output txt file (default: "GC_out.txt")

required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        fasta file
```

### Calculate GC and sliding window

```
$ bit-calc-gc-sliding-window -h
usage: bit-calc-gc-sliding-window [-h] [-i INPUT_FASTA]
                                             [-o OUTPUT_FILE] [-w WINDOW]
                                             [-s STEP]

This script is for nucleotide multifastas and will return a tab-delimited file
with 4 columns: header, sequence length, gc of whole sequence, and gc of each
window of the specified window size (-w) for each step of the specified step
size (-s).

optional arguments:
  -h, --help            show this help message and exit
  -w WINDOW, --window_size WINDOW
                        Desired size of sliding window (default: 100)
  -s STEP, --step_size STEP
                        Desired size of steps between each window (default: 1)

required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        fasta file
  -o OUTPUT_FILE, --output_txt_file OUTPUT_FILE
                        Name of output txt file
```

### Extract sequences by coordinates

```
$ bit-extract-seqs-by-coords -h
usage: bit-extract-seqs-by-coords [-h] [-i INPUT_FASTA] [-b BED_FILE]
                                       [-o OUTPUT_FASTA]

This script takes a multifasta file and tab-delimited file specifying which
contigs and coordinates are wanted and returns a multifasta of the chopped out
sequences. NOTE: It requires the python package "pybedtools".

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FASTA, --output_fasta OUTPUT_FASTA
                        Name of output fasta file (default:
                        "Extracted_seqs.fa")

required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Starting fasta file
  -b BED_FILE, --bed_file BED_FILE
                        Bed file of desired contigs and coordinates (3 columns
                        - contig, start, end - no header, 0-based counting
```

### Parse fasta by headers

```
$ bit-parse-fasta-by-headers -h
usage: parse-fasta-by-headers [-h] -i INPUT_FASTA -w HEADERS [-o OUTPUT_FASTA]
                              [--inverse]

This script is for parsing a fasta file by pulling out sequences with the
desired headers. If you want all sequences EXCEPT the ones with the headers
you are providing, add the flag "--inverse".

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FASTA, --output_fasta OUTPUT_FASTA
                        Output fasta file default: "Wanted.fa"
  --inverse             Add this flag to pull out all sequences with headers
                        NOT in the provided header file.

required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Original fasta file
  -w HEADERS, --sequence_headers HEADERS
                        Single-column file with sequence headers
```

### Simplify sequence headers

```
$ bit-simplify-fasta-headers -h
usage: bit-simplify-fasta-headers [-h] -i INPUT_FASTA [-w WANTED_NAME]
                                  [-o OUTPUT_PREFIX]

This script will rename all sequences of a multifasta with the same name with
an appended number to keep them unique.

optional arguments:
  -h, --help            show this help message and exit
  -w WANTED_NAME, --desired_name WANTED_NAME
                        Name to give seqs (default: "Seq"
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix added to output fasta (default: "Renamed").

required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Starting fasta file
```

### Reorder fasta by headers

```
$ bit-reorder-fasta -h
usage: This script takes a multifasta file and reorders the sequences according to the headers provided.
       [-h] -i INPUT_FASTA -w ORDERED_HEADERS [-o OUTPUT_PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix added to output fasta (default: "Reordered").

required arguments:
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        Original fasta file
  -w ORDERED_HEADERS, --wanted_sequence_order ORDERED_HEADERS
                        Single-column file with headers in desired order
```

### Calculator

```
$ bit-calc -h
This script echoes the input into the bash program `bc` and returns the answer,
saving very little time compared to how much I waste.

Usage:
	 bit-calc "(5+5)/2"
```
