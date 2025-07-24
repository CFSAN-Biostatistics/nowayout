<p align="center">
    <img src="../assets/nowayout-icon.png" width="20%" height="20%" />
</p>

---

`nowayout` is a **super-fast** automated software pipeline for taxonomic classification of Eukaryotic mitochondrial reads. It uses a custom database to first identify mitochondrial reads and performs read classification on those identified reads.

---

<!-- TOC -->

- [Minimum Requirements](#minimum-requirements)
- [HFP GalaxyTrakr](#hfp-galaxytrakr)
- [Usage and Examples](#usage-and-examples)
  - [Databases](#databases)
  - [Input](#input)
  - [Output](#output)
  - [Preset filters](#preset-filters)
  - [Computational resources](#computational-resources)
  - [Runtime profiles](#runtime-profiles)
  - [your_institution.config](#your_institutionconfig)
- [Test Run](#test-run)
- [nowayout CLI Help](#nowayout-cli-help)

<!-- /TOC -->

\
&nbsp;

## Minimum Requirements

1. [Nextflow version 25.04.6](https://github.com/nextflow-io/nextflow/releases/download/v25.04.6/nextflow).
    - Make the `nextflow` binary executable (`chmod 755 nextflow`) and also make sure that it is made available in your `$PATH`.
    - If your existing `JAVA` install does not support the newest **Nextflow** version, you can try **Amazon**'s `JAVA` (OpenJDK):  [Corretto](https://docs.aws.amazon.com/corretto/latest/corretto-21-ug/downloads-list.html).
2. Either of `micromamba` (version `1.5.9`) or `docker` or `singularity` installed and made available in your `$PATH`.
    - Running the workflow via `micromamba` software provisioning is **preferred** as it does not require any `sudo` or `admin` privileges or any other configurations with respect to the various container providers.
    - To install `micromamba` for your system type, please follow these [installation steps](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html#linux-and-macos) and make sure that the `micromamba` binary is made available in your `$PATH`.
    - Just the `curl` step is sufficient to download the binary as far as running the workflows are concerned.
    - Once you have finished the installation, **it is important that you downgrade `micromamba` to version `1.5.9`**.
    - First check, if your version is other than `1.5.9` and if not, do the downgrade.

        ```bash
        micromamba --version
        micromamba self-update --version 1.5.9 -c conda-forge
        ```

3. Minimum of 10 CPU cores and about 60 GBs for main workflow steps. More memory may be required if your **FASTQ** files are big.

\
&nbsp;

## HFP GalaxyTrakr

The `nowayout` pipeline **will** be made available for use on the newest version of [Galaxy instance supported by HFP, FDA](https://galaxytrakr.org/) (`version >= 24.x`). Please check this space for announcements in this regard.

Please note that the pipeline on [HFP GalaxyTrakr](https://galaxytrakr.org) in most cases may be a version older than the one on **GitHub** due to testing prioritization.

\
&nbsp;

## Usage and Examples

Clone or download this repository and then call `cpipes`.

```bash
cpipes --pipeline nowayout [options]
```

Alternatively, you can use `nextflow` to directly pull and run the pipeline.

```bash
nextflow pull CFSAN-Biostatistics/nowayout
nextflow list
nextflow info CFSAN-Biostatistics/nowayout
nextflow run CFSAN-Biostatistics/nowayout --pipeline nowayout --help
```

\
&nbsp;

### Databases

---

The successful run of the workflow requires proper setup of the custom database files:

- `nowayout_dbs`: [Download](https://cfsan-pub-xfer.s3.amazonaws.com/Kranti.Konganti/nowayout/nowayout_dbs.tar.bz2) (~ 22 GB).

Once you have downloaded the databases, uncompress and set the **UNIX symbolic link** to the database folders in [assets](../assets/) folder as follows:

```bash
mkdir assets/dbfiles
cd assets/dbfiles
ln -s /path/to/nowayout_dbs/kma kma
ln -s /path/to/nowayout_dbs/reference reference
ln -s /path/to/nowayout_dbs/taxonomy taxonomy
```

That's it!

\
&nbsp;

### Input

---

The input to the workflow is a folder containing compressed (`.gz`) FASTQ files of long reads or short reads. Please note that the sample grouping happens automatically by the file name of the FASTQ file. If for example, a single sample is sequenced across multiple sequencing lanes, you can choose to group those FASTQ files into one sample by using the `--fq_filename_delim` and `--fq_filename_delim_idx` options. By default, `--fq_filename_delim` is set to `_` (underscore) and `--fq_filename_delim_idx` is set to 1.

For example, if the directory contains FASTQ files as shown below:

- KB-01_apple_L001_R1.fastq.gz
- KB-01_apple_L001_R2.fastq.gz
- KB-01_apple_L002_R1.fastq.gz
- KB-01_apple_L002_R2.fastq.gz
- KB-02_mango_L001_R1.fastq.gz
- KB-02_mango_L001_R2.fastq.gz
- KB-02_mango_L002_R1.fastq.gz
- KB-02_mango_L002_R2.fastq.gz

Then, to create 2 sample groups, `apple` and `mango`, we split the file name by the delimitor (underscore in the case, which is default) and group by the first 2 words (`--fq_filename_delim_idx 2`).

This goes without saying that all the FASTQ files should have uniform naming patterns so that `--fq_filename_delim` and `--fq_filename_delim_idx` options do not have any adverse effect in collecting and creating a sample metadata sheet.

\
&nbsp;

### Output

---

All the outputs for each step are stored inside the folder mentioned with the `--output` option. A `multiqc_report.html` file inside the `nowayout-multiqc` folder can be opened in any browser on your local workstation which contains a consolidated brief report.

Please note that the percentage relative abundances seen are relative to the total number of mitochondrial reads and not the total number of reads per sample.

\
&nbsp;

### Preset filters

---

There are three preset threshold filters that are available with the pipeline: `--nowo_thresholds strict`, `--nowo_thresholds mild` and `--nowo_thresholds relax`. Use these options for exploration of results via multiple runs on the same input dataset. The default is `strict` thresholds.

\
&nbsp;

### Computational resources

---

The workflows `nowayout` require at least a minimum of 10 CPU cores and 60 GBs of memory to successfully finish the workflow.

\
&nbsp;

### Runtime profiles

---

You can use different run time profiles that suit your specific compute environments i.e., you can run the workflow locally on your machine or in a grid computing infrastructure.

\
&nbsp;

Example:

```bash
cd /data/scratch/$USER
mkdir nf-cpipes
cd nf-cpipes
cpipes \
    --pipeline nowayout \
    --input /path/to/fastq_pass_dir \
    --output /path/to/where/output/should/go \
    -profile your_institution
```

The above command would run the pipeline and store the output at the location per the `--output` flag and the **NEXTFLOW** reports are always stored in the current working directory from where `cpipes` is run. For example, for the above command, a directory called `CPIPES-nowayout` would hold all the **NEXTFLOW** related logs, reports and trace files.

\
&nbsp;

### `your_institution.config`

---

In the above example, we can see that we have mentioned the run time profile as `your_institution`. For this to work, add the following lines at the end of [`computeinfra.config`](../conf/computeinfra.config) file which should be located inside the `conf` folder. For example, if your institution uses **SGE** or **UNIVA** for grid computing instead of **SLURM** and has a job queue named `normal.q`, then add these lines:

\
&nbsp;

```groovy
your_institution {
    process.executor = 'sge'
    process.queue = 'normal.q'
    singularity.enabled = false
    singularity.autoMounts = true
    docker.enabled = false
    params.enable_conda = true
    conda.enabled = true
    conda.useMicromamba = true
    params.enable_module = false
}
```

In the above example, by default, all the software provisioning choices are disabled except `conda`. You can also choose to remove the `process.queue` line altogether and the `nowayout` workflow will request the appropriate memory and number of CPU cores automatically, which ranges from 1 CPU, 1 GB and 1 hour for job completion up to 10 CPU cores, 1 TB and 120 hours for job completion.

\
&nbsp;

### Cloud computing

---

You can run the workflow in the cloud (works only with proper set up of AWS resources). Add new run time profiles with required parameters per [Nextflow docs](https://www.nextflow.io/docs/latest/executor.html):

\
&nbsp;

Example:

```groovy
my_aws_batch {
    executor = 'awsbatch'
    queue = 'my-batch-queue'
    aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'
    aws.batch.region = 'us-east-1'
    singularity.enabled = false
    singularity.autoMounts = true
    docker.enabled = true
    params.conda_enabled = false
    params.enable_module = false
}
```

\
&nbsp;

## Test Run

After you make sure that you have all the [minimum requirements](#minimum-requirements) to run the workflow, you can try the `nowayout` on some datasets.

- Download input reads [from S3](https://cfsan-pub-xfer.s3.amazonaws.com/Kranti.Konganti/nowayout/nowayout_test_reads.tar.bz2) (~ 8 GB).
  - This dataset was part of the research for detecting and identifying insects or insect fragments in food, an essential component of food safety and regulatory monitoring. Insects such as **_Plodia interpunctella_** (Indian meal moth) and _**Tribolium castaneum**_ (red flour beetle) were intentionally spiked into wheat flour at varying concentrations to create benchmark samples. These serve as reference materials to test and validate molecular detection workflows.
- Download pre-formatted  databases (**MANDATORY**) [from S3](https://cfsan-pub-xfer.s3.amazonaws.com/Kranti.Konganti/nowayout/nowayout_dbs.tar.bz2) (~ 22 GB).
- After successful download, untar and add **symbolic links** in [assets](../assets) folder as described in the [Databases](#databases) section.
- It is always a best practice to use absolute UNIX paths and real destinations of symbolic links during pipeline execution. For example, find out the real path(s) of your absolute UNIX path(s) and use that for the `--input` and `--output` options of the pipeline.

  ```bash
  realpath /hpc/scratch/user/input/srr
  ```

- Now run the workflow by ignoring quality values since these are simulated base qualities:

    ```bash
    cpipes \
        --pipeline nowayout \
        --input /path/to/nowayout_test_reads \
        --output /path/to/nowayout_test_output \
        --fq_single_end true \
        -profile stdkondagac \
        -resume
    ```

- After succesful run of the workflow, your **MultiQC** report should look something like [this](https://cfsan-pub-xfer.s3.us-east-1.amazonaws.com/Kranti.Konganti/nowayout/CPIPES-Report_multiqc_report.html).

- `nowayout` also automatically generates [Krona](https://github.com/marbl/Krona) charts. The **Krona** chart for the above test run should look something like [this](https://cfsan-pub-xfer.s3.us-east-1.amazonaws.com/Kranti.Konganti/nowayout/CPIPES_nowayout_krona.html)

Please note that the run time profile `stdkondagac` will run jobs locally using `micromamba` for software provisioning. The first time you run the command, a new folder called `kondagac_cache` will be created and subsequent runs should use this `conda` cache.

\
&nbsp;

## `nowayout` CLI Help

```text
cpipes --pipeline nowayout --help

 N E X T F L O W   ~  version 24.10.4

Launching `/home/user/nowayout/cpipes` [sleepy_pauling] DSL2 - revision: 55d6f63710

================================================================================
             (o)                  
  ___  _ __   _  _ __    ___  ___ 
 / __|| '_ \ | || '_ \  / _ \/ __|
| (__ | |_) || || |_) ||  __/\__ \
 \___|| .__/ |_|| .__/  \___||___/
      | |       | |               
      |_|       |_|
--------------------------------------------------------------------------------
A collection of modular pipelines at CFSAN, FDA.
--------------------------------------------------------------------------------
Name                            : CPIPES
Author                          : Kranti.Konganti@fda.hhs.gov
Version                         : 0.8.0
Center                          : CFSAN, FDA.
================================================================================

Workflow                        : nowayout

Author                          : Kranti Konganti

Version                         : 0.5.0


Usage                           : cpipes --pipeline nowayout [options]


Required                        : 

--input                         : Absolute path to directory containing FASTQ
                                  files. The directory should contain only
                                  FASTQ files as all the files within the
                                  mentioned directory will be read. Ex: --
                                  input /path/to/fastq_pass

--output                        : Absolute path to directory where all the
                                  pipeline outputs should be stored. Ex: --
                                  output /path/to/output

Other options                   : 

--metadata                      : Absolute path to metadata CSV file
                                  containing five mandatory columns: sample,
                                  fq1,fq2,strandedness,single_end. The fq1
                                  and fq2 columns contain absolute paths to
                                  the FASTQ files. This option can be used in
                                  place of --input option. This is rare. Ex
                                  : --metadata samplesheet.csv

--fq_suffix                     : The suffix of FASTQ files (Unpaired reads
                                  or R1 reads or Long reads) if an input
                                  directory is mentioned via --input option.
                                  Default: _R1_001.fastq.gz

--fq2_suffix                    : The suffix of FASTQ files (Paired-end reads
                                  or R2 reads) if an input directory is
                                  mentioned via --input option. Default:
                                  _R2_001.fastq.gz

--fq_filter_by_len              : Remove FASTQ reads that are less than this
                                  many bases. Default: 0

--fq_strandedness               : The strandedness of the sequencing run.
                                  This is mostly needed if your sequencing
                                  run is RNA-SEQ. For most of the other runs
                                  , it is probably safe to use unstranded for
                                  the option. Default: unstranded

--fq_single_end                 : SINGLE-END information will be auto-
                                  detected but this option forces PAIRED-END
                                  FASTQ files to be treated as SINGLE-END so
                                  only read 1 information is included in auto
                                  -generated samplesheet. Default: false

--fq_filename_delim             : Delimiter by which the file name is split
                                  to obtain sample name. Default: _

--fq_filename_delim_idx         : After splitting FASTQ file name by using
                                  the --fq_filename_delim option, all
                                  elements before this index (1-based) will
                                  be joined to create final sample name.
                                  Default: 1

--fastp_run                     : Run fastp tool. Default: true

--fastp_failed_out              : Specify whether to store reads that cannot
                                  pass the filters. Default: false

--fastp_merged_out              : Specify whether to store merged output or
                                  not. Default: false

--fastp_overlapped_out          : For each read pair, output the overlapped
                                  region if it has no mismatched base.
                                  Default: false

--fastp_6                       : Indicate that the input is using phred64
                                  scoring (it'll be converted to phred33, so
                                  the output will still be phred33). Default
                                  : false

--fastp_reads_to_process        : Specify how many reads/pairs are to be
                                  processed. Default value 0 means process
                                  all reads. Default: 0

--fastp_fix_mgi_id              : The MGI FASTQ ID format is not compatible
                                  with many BAM operation tools, enable this
                                  option to fix it. Default: false

--fastp_A                       : Disable adapter trimming. On by default.
                                  Default: false

--fastp_adapter_fasta           : Specify a FASTA file to trim both read1 and
                                  read2 (if PE) by all the sequences in this
                                  FASTA file. Default: false

--fastp_f                       : Trim how many bases in front of read1.
                                  Default: 0

--fastp_t                       : Trim how many bases at the end of read1.
                                  Default: 0

--fastp_b                       : Max length of read1 after trimming. Default
                                  : 0

--fastp_F                       : Trim how many bases in front of read2.
                                  Default: 0

--fastp_T                       : Trim how many bases at the end of read2.
                                  Default: 0

--fastp_B                       : Max length of read2 after trimming. Default
                                  : 0

--fastp_dedup                   : Enable deduplication to drop the duplicated
                                  reads/pairs. Default: true

--fastp_dup_calc_accuracy       : Accuracy level to calculate duplication (1~
                                  6), higher level uses more memory (1G, 2G,
                                  4G, 8G, 16G, 24G). Default 1 for no-dedup
                                  mode, and 3 for dedup mode. Default: 6

--fastp_poly_g_min_len          : The minimum length to detect polyG in the
                                  read tail. Default: 10

--fastp_G                       : Disable polyG tail trimming. Default: true

--fastp_x                       : Enable polyX trimming in 3' ends. Default:
                                  false

--fastp_poly_x_min_len          : The minimum length to detect polyX in the
                                  read tail. Default: 10

--fastp_cut_front               : Move a sliding window from front (5') to
                                  tail, drop the bases in the window if its
                                  mean quality < threshold, stop otherwise.
                                  Default: true

--fastp_cut_tail                : Move a sliding window from tail (3') to
                                  front, drop the bases in the window if its
                                  mean quality < threshold, stop otherwise.
                                  Default: false

--fastp_cut_right               : Move a sliding window from tail, drop the
                                  bases in the window and the right part if
                                  its mean quality < threshold, and then stop
                                  . Default: true

--fastp_W                       : Sliding window size shared by --
                                  fastp_cut_front, --fastp_cut_tail and --
                                  fastp_cut_right. Default: 20

--fastp_M                       : The mean quality requirement shared by --
                                  fastp_cut_front, --fastp_cut_tail and --
                                  fastp_cut_right. Default: 30

--fastp_q                       : The quality value below which a base should
                                  is not qualified. Default: 30

--fastp_u                       : What percent of bases are allowed to be
                                  unqualified. Default: 40

--fastp_n                       : How many N's can a read have. Default: 5

--fastp_e                       : If the full reads' average quality is below
                                  this value, then it is discarded. Default
                                  : 0

--fastp_l                       : Reads shorter than this length will be
                                  discarded. Default: 35

--fastp_max_len                 : Reads longer than this length will be
                                  discarded. Default: 0

--fastp_y                       : Enable low complexity filter. The
                                  complexity is defined as the percentage of
                                  bases that are different from its next base
                                  (base[i] != base[i+1]). Default: true

--fastp_Y                       : The threshold for low complexity filter (0~
                                  100). Ex: A value of 30 means 30%
                                  complexity is required. Default: 30

--fastp_U                       : Enable Unique Molecular Identifier (UMI)
                                  pre-processing. Default: false

--fastp_umi_loc                 : Specify the location of UMI, can be one of
                                  index1/index2/read1/read2/per_index/
                                  per_read. Default: false

--fastp_umi_len                 : If the UMI is in read1 or read2, its length
                                  should be provided. Default: false

--fastp_umi_prefix              : If specified, an underline will be used to
                                  connect prefix and UMI (i.e. prefix=UMI,
                                  UMI=AATTCG, final=UMI_AATTCG). Default:
                                  false

--fastp_umi_skip                : If the UMI is in read1 or read2, fastp can
                                  skip several bases following the UMI.
                                  Default: false

--fastp_p                       : Enable overrepresented sequence analysis.
                                  Default: true

--fastp_P                       : One in this many number of reads will be
                                  computed for overrepresentation analysis (1
                                  ~10000), smaller is slower. Default: 20

--kmaalign_run                  : Run kma tool. Default: true

--kmaalign_int                  : Input file has interleaved reads.  Default
                                  : false

--kmaalign_ef                   : Output additional features. Default: false

--kmaalign_vcf                  : Output vcf file. 2 to apply FT. Default:
                                  false

--kmaalign_sam                  : Output SAM, 4/2096 for mapped/aligned.
                                  Default: false

--kmaalign_nc                   : No consensus file. Default: true

--kmaalign_na                   : No aln file. Default: true

--kmaalign_nf                   : No frag file. Default: false

--kmaalign_a                    : Output all template mappings. Default:
                                  false

--kmaalign_and                  : Use both -mrs and p-value on consensus.
                                  Default: true

--kmaalign_oa                   : Use neither -mrs or p-value on consensus.
                                  Default: false

--kmaalign_bc                   : Minimum support to call bases. Default:
                                  false

--kmaalign_bcNano               : Altered indel calling for ONT data. Default
                                  : false

--kmaalign_bcd                  : Minimum depth to call bases. Default: false

--kmaalign_bcg                  : Maintain insignificant gaps. Default: false

--kmaalign_ID                   : Minimum consensus ID. Default: 85.0

--kmaalign_md                   : Minimum depth. Default: false

--kmaalign_dense                : Skip insertion in consensus. Default: false

--kmaalign_ref_fsa              : Use Ns on indels. Default: false

--kmaalign_Mt1                  : Map everything to one template. Default:
                                  false

--kmaalign_1t1                  : Map one query to one template. Default:
                                  false

--kmaalign_mrs                  : Minimum relative alignment score. Default:
                                  0.99

--kmaalign_mrc                  : Minimum query coverage. Default: 0.99

--kmaalign_mp                   : Minimum phred score of trailing and leading
                                  bases. Default: 30

--kmaalign_mq                   : Set the minimum mapping quality. Default:
                                  false

--kmaalign_eq                   : Minimum average quality score. Default: 30

--kmaalign_5p                   : Trim 5 prime by this many bases. Default:
                                  false

--kmaalign_3p                   : Trim 3 prime by this many bases Default:
                                  false

--kmaalign_apm                  : Sets both -pm and -fpm Default: false

--kmaalign_cge                  : Set CGE penalties and rewards Default:
                                  false

--seqkit_grep_run               : Run the seqkit `grep` tool. Default: true

--seqkit_grep_n                 : Match by full name instead of just ID.
                                  Default: undefined

--seqkit_grep_s                 : Search subseq on seq, both positive and
                                  negative strand are searched, and mismatch
                                  allowed using flag --seqkit_grep_m. Default
                                  : undefined

--seqkit_grep_c                 : Input is circular genome Default: undefined

--seqkit_grep_C                 : Just print a count of matching records.
                                  With the --seqkit_grep_v flag, count non-
                                  matching records. Default: undefined

--seqkit_grep_i                 : Ignore case while using seqkit grep.
                                  Default: undefined

--seqkit_grep_v                 : Invert the match i.e. select non-matching
                                  records. Default: undefined

--seqkit_grep_m                 : Maximum mismatches when matching by
                                  sequence. Default: undefined

--seqkit_grep_r                 : Input patters are regular expressions.
                                  Default: undefined

--salmonidx_run                 : Run `salmon index` tool. Default: true

--salmonidx_k                   : The size of k-mers that should be used for
                                  the  quasi index. Default: false

--salmonidx_gencode             : This flag will expect the input transcript
                                  FASTA to be in GENCODE format, and will
                                  split the transcript name at the first `|`
                                  character. These reduced names will be used
                                  in the output and when looking for these
                                  transcripts in a gene to transcript GTF.
                                  Default: false

--salmonidx_features            : This flag will expect the input reference
                                  to be in the tsv file format, and will
                                  split the feature name at the first `tab`
                                  character. These reduced names will be used
                                  in the output and when looking for the
                                  sequence of the features. GTF. Default:
                                  false

--salmonidx_keepDuplicates      : This flag will disable the default indexing
                                  behavior of discarding sequence-identical
                                  duplicate transcripts. If this flag is
                                  passed then duplicate transcripts that
                                  appear in the input will be retained and
                                  quantified separately. Default: true

--salmonidx_keepFixedFasta      : Retain the fixed fasta file (without short
                                  transcripts and duplicates, clipped, etc.)
                                  generated during indexing. Default: false

--salmonidx_filterSize          : The size of the Bloom filter that will be
                                  used by TwoPaCo during indexing. The filter
                                  will be of size 2^{filterSize}. A value of
                                  -1 means that the filter size will be
                                  automatically set based on the number of
                                  distinct k-mers in the input, as estimated
                                  by nthll. Default: false

--salmonidx_sparse              : Build the index using a sparse sampling of
                                  k-mer positions This will require less
                                  memory (especially during quantification),
                                  but will take longer to constructand can
                                  slow down mapping / alignment. Default:
                                  false

--salmonidx_n                   : Do not clip poly-A tails from the ends of
                                  target sequences. Default: true

--sourmashsketch_run            : Run `sourmash sketch dna` tool. Default:
                                  true

--sourmashsketch_mode           : Select which type of signatures to be
                                  created: dna, protein, fromfile or
                                  translate. Default: dna

--sourmashsketch_p              : Signature parameters to use. Default: '
                                  abund,scaled=100,k=71

--sourmashsketch_file           : <path>  A text file containing a list of
                                  sequence files to load. Default: false

--sourmashsketch_f              : Recompute signatures even if the file
                                  exists. Default: false

--sourmashsketch_name           : Name the signature generated from each file
                                  after the first record in the file.
                                  Default: false

--sourmashsketch_randomize      : Shuffle the list of input files randomly.
                                  Default: false

--sourmashgather_run            : Run `sourmash gather` tool. Default: true

--sourmashgather_n              : Number of results to report. By default,
                                  will terminate at --sourmashgather_thr_bp
                                  value. Default: false

--sourmashgather_thr_bp         : Reporting threshold (in bp) for estimated
                                  overlap with remaining query. Default: 100

--sourmashgather_ani_ci         : Output confidence intervals for ANI
                                  estimates. Default: true

--sourmashgather_k              : The k-mer size to select. Default: 71

--sourmashgather_dna            : Choose DNA signature. Default: true

--sourmashgather_rna            : Choose RNA signature. Default: false

--sourmashgather_nuc            : Choose Nucleotide signature. Default: false

--sourmashgather_scaled         : Scaled value should be between 100 and 1e6
                                  . Default: false

--sourmashgather_inc_pat        : Search only signatures that match this
                                  pattern in name, filename, or md5. Default
                                  : false

--sourmashgather_exc_pat        : Search only signatures that do not match
                                  this pattern in name, filename, or md5.
                                  Default: false

--sfhpy_run                     : Run the sourmash_filter_hits.py script.
                                  Default: true

--sfhpy_fcn                     : Column name by which filtering of rows
                                  should be applied. Default: f_match

--sfhpy_fcv                     : Remove genomes whose match with the query
                                  FASTQ is less than this much. Default: 0.8

--sfhpy_gt                      : Apply greather than or equal to condition
                                  on numeric values of --sfhpy_fcn column.
                                  Default: true

--sfhpy_lt                      : Apply less than or equal to condition on
                                  numeric values of --sfhpy_fcn column.
                                  Default: false

--sfhpy_all                     : Instead of just the column value, print
                                  entire row. Default: true

--gsalkronapy_run               : Run the `gen_salmon_tph_and_krona_tsv.py`
                                  script. Default: true

--gsalkronapy_sf                : Set the scaling factor by which TPM values
                                  are scaled down. Default: 10000

--gsalkronapy_smres_suffix      : Find the `sourmash gather` result files
                                  ending in this suffix. Default: false

--gsalkronapy_failed_suffix     : Find the sample names which failed
                                  classification stored inside the files
                                  ending in this suffix. Default: false

--gsalkronapy_num_lin_cols      : Number of columns expected in the lineages
                                  CSV file.  Default: false

--gsalkronapy_lin_regex         : Number of columns expected in the lineages
                                  CSV file.  Default: false

--krona_ktIT_run                : Run the ktImportText (ktIT) from krona.
                                  Default: true

--krona_ktIT_n                  : Name of the highest level. Default: all

--krona_ktIT_q                  : Input file(s) do not have a field for
                                  quantity. Default: false

--krona_ktIT_c                  : Combine data from each file, rather than
                                  creating separate datasets within the chart
                                  . Default: false

Help options                    : 

--help                          : Display this message.
```
