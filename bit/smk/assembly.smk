import os
from pathlib import Path

reads_dict = config['reads_json']
output_dir = config['output_dir']
base_output_dir = config['base_output_dir']
log_files_dir = config['log_files_dir']
base_log_files_dir = config['base_log_files_dir']
reads_dir = config['reads_dir']
assembler = config['assembler']
threads = config['threads']
skip_fastp = config['skip_fastp']
skip_bbnorm = config['skip_bbnorm']
memory = config['memory']
min_contig_len = config['min_contig_len']
isolate = config['isolate']

print(reads_dict)

