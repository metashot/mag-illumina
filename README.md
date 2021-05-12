# mag-illumina

metashot/mag-illumina is a workflow for the assembly and binning of Illumina
sequences from metagenomic samples.

- [MetaShot Home](https://metashot.github.io/)

## Main features

- Input: single-end, paired-end (also interleaved) Illumina sequences (gzip and
  bzip2 compressed FASTQ also supported);
- Histogram text files (for each input sample) of base frequency, quality
  scores, GC content, average quality and length are generated from input reads
  and clean reads using
  [bbduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/);
- Adapter trimming, contaminant filtering and quality filtering/trimming and
  length filtering using bbduk;
- Assembly with [Spades](https://cab.spbu.ru/software/spades/) or
  [Megahit](https://github.com/voutcn/megahit);
- Assembly statistics using
  [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/statistics-guide/);
- Binning with
  [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/).
- Optionally, assemble plasmids with metaplasmidSPAdes and verify them using
  ViralVerify.

## Quick start

1. Install Docker (or Singulariry) and Nextflow (see
   [Dependences](https://metashot.github.io/#dependencies));
1. Start running the analysis:
   
  ```bash
  nextflow run metashot/mag-illumina \
    --reads '*_R{1,2}.fastq.gz' \
    --outdir results
  ```

## Parameters
See the file [`nextflow.config`](nextflow.config) for the complete list of
parameters.

## Output
The files and directories listed below will be created in the `results` directory
after the pipeline has finished.

### Main outputs

- `scaffolds`: scaffolds for each input sample;
- `bins`: genome bins;
- `unbinned`: unbinned contigs;
- `stats_scaffolds.tsv`: scaffold statistics;
- `verified_plasmids`: verified plasmids (if `--run_metaplasmidspades` is set).

### Secondary outputs

- `raw_reads_stats`: base frequency, quality scores, gc content, average
  quality and length for each input sample;
- `clean_reads_stats`: same as above, but for the reads after the quality
  control;
- `clean_reads`: clean reads (if `--save_clean` is set);
- `qc`: adapter trimming and contaminant filtering statistics;
- `metaspades`, `metaplasmidspades` and `megahit`: complete assembler output for
  each sample (if `--save_assembler_output` is set);
- `scaffolds_plasmids`: candidate plasmids (if `--run_metaplasmidspades` is set);
- `viralverify`: viralVerify output (if `--run_metaplasmidspades` is set);
- `metabat2`: metabat2 log and the depth of coverage for each assembly.

## System requirements
Please refer to [System
requirements](https://metashot.github.io/#system-requirements) for the complete
list of system requirements options.
