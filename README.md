# de-novo-assembly
Automated de novo assembly of selected genes using BBMap and Trinity
The SRA file is downloaded with fastq-dump program available from SRA-Toolkit (https://hpc.nih.gov/apps/sratoolkit.html). The reads aligned to sequence of interest are selected with BBMap program (by Bushnell B.  https://sourceforge.net/projects/bbmap/) in ‘vslow’ and ‘local’ modes and “maxindel” set to 100. Next, the paired reads in the fastq file exported by BBMap are separated into two fastq files with bbsplitpairs scripts from BBMap, which are then assembled de novo by Trinity (https://github.com/trinityrnaseq/trinityrnaseq) three separate times as follow: (1) --KMER_SIZE 32, (2) stringent setting, which includes “--min_kmer_cov 4 –min_glue 4 –min_iso_ratio 0.2 –glue_factor 0.2 –jaccard_clip”, (3) both –KMER_SIZE 32 and stringent setting. If there are more than 10,000 reads in each fastq file, the first 5000 reads extracted with seqtk program (https://github.com/lh3/seqtk) are assembled in two more Trinity runs with –KMER_SIZE 32 with or without the stringent setting. The read coverages of starting bases are than obtained for assemblies that covered at least 90% of the reference sequences with alignment scores greater than 350 using  BBMap under “perfectmode” and “startcov=t”. We automated the above process with python scripts, which can be executed in Windows Subsystem for Linux from a shell script file, which can be supplied with multiple SRA IDs for high-throughput assembly.
Instructions for how to run these scripts are included in "run_get_coverage.sh".
