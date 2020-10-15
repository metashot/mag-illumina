# metashot/mag-illumina Nextflow

metashot/mag-illumina is a [Nextflow](https://www.nextflow.io/) pipeline for
the assembly and binning of Illumina sequences from metagenomic samples.

Main features:

- Input: single-end, paired-end (also interleaved) Illumina sequences (gzip
  and bzip2 compressed FASTQ also supported);
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

## Quick start

1. Install [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/);
1. Start running the analysis:
   
  ```bash
  nextflow run metashot/mag-illumina
    --input '*_R{1,2}.fastq.gz' \
    --outdir results
  ```

See the file [`nextflow.config`](nextflow.config) for the complete list of parameters.

### Output
Several directories will be created in the `results` folder:

- `raw_reads_stats`: base frequency, quality scores, gc content, average
  quality and length for each input sample;
- `clean_reads_stats`: same as above, but for the reads after the quality
  control;
- `qc`: statistaics about the adapter trimming and the contaminant filtering;
- `spades`: Spades log for each input sample;
- `scaffolds`: scaffolds for each input samples;
- `assembly_stats`: scaffold statistics for each each assembly;
- `metabat2`: metabat2 log and the depth of coverage for each assembly;
- `bins`: metabat2 bins for each assembly;
- `unbinned`: unbinned contigs for each assembly.

## System requirements
Each step in the pipeline has a default set of requirements for number of CPUs,
memory and time. For some of the steps in the pipeline, if the job exits with an
error it will automatically resubmit with higher requests (see
[`process.config`](process.config)).

You can customize the compute resources that the pipeline requests by either:
- setting the global parameters `--max_cpus`, `--max_memory` and
  `--max_time`, or
- creating a [custom config
  file](https://www.nextflow.io/docs/latest/config.html#configuration-file)
  (`-c` or `-C` parameters), or
- modifying the [`process.config`](process.config) file.

## Reproducibility
We recommend to specify a pipeline version when running the pipeline on your
data with the `-r` parameter:

```bash
  nextflow run metashot/kraken2 -r 1.0.0
    ...
```

Moreover, this workflow uses the docker images available at
https://hub.docker.com/u/metashot/ for reproducibility. You can check the
version of the software used in the workflow by opening the file
[`process.config`](process.config). For example `container =
metashot/kraken2:2.0.9-beta-6` means that the version of kraken2 is the
`2.0.9-beta` (the last number, 6, is the metashot release of this container).

## Singularity
If you want to use [Singularity](https://singularity.lbl.gov/) instead of Docker,
comment the Docker lines in [`nextflow.config`](nextflow.config) and add the following:

```nextflow
singularity.enabled = true
singularity.autoMounts = true
```

## Credits
This workflow is maintained Davide Albanese and Claudio Donati at the [FEM's
Unit of Computational
Biology](https://www.fmach.it/eng/CRI/general-info/organisation/Chief-scientific-office/Computational-biology).
