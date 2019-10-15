# nanopore_reads_binner

Bin MinION reads by (cumulative) intervals.

## Usage
```
usage: nanopore_reads_binner.py [-h] -f /basecalled/folder/ -o /qc/ [-t 48] -i
                                1h [-p my_output_name]

Bin MinION reads by (cumulative) intervals

optional arguments:
  -h, --help            show this help message and exit
  -f /basecalled/folder/, --fastq /basecalled/folder/
                        Input folder with fastq file(s), gzipped or not
  -o /qc/, --output /qc/
                        Output folder
  -t 48, --threads 48   Number of CPUDefault: 48
  -i 1h, --interval 1h  Time interval to create the bins.Use integer values
                        (no decimal)Can use 'h', 'm' or 's'Bins are
                        cumulative. If interval is set to at 1h, second bin
                        willalso contains the reads from the first
                        bin.Default: 1h
  -p my_output_name, --prefix my_output_name
                        Output file prefixIf using "my_sample" as prefix,
                        files will be named "my_sample_1h.fastq.gz",
                        etc.Default: interval

  
```
