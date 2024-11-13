- **run_kraken_braken.sh:** bash script to extract unaligned reads, realign to CHM13, then re-extract unaligned reads and trim. Next, classify with kraken and bracken. Requires 3 inputs in order: the sequencing type (WGS,RNA,16S), the name of the sample for outputting, and lastly the location of the input bam, OR if handling 16S where input is expected to be unaligned, two paired fastq files.
