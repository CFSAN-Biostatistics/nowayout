# NextFlow DSL2 Module

```bash
SEQKIT_GREP
```

## Description

Run `seqkit grep` command on reads in FASTQ format. Produces a filtered FASTQ file as per the filter strategy in the supplied input file.

\
&nbsp;

### `input:`

___

Type: `tuple`

Takes in the following tuple of metadata (`meta`) and a list of reads of type `path` (`reads`) per sample (`id:`).

Ex:

```groovy
[
    [ id: 'FAL00870',
       strandedness: 'unstranded',
       single_end: true,
       centrifuge_x: '/hpc/db/centrifuge/2022-04-12/ab'
    ],
    '/hpc/scratch/test/FAL000870/f1.merged.fq.gz'
]
```

\
&nbsp;

#### `meta`

Type: Groovy Map

A Groovy Map containing the metadata about the FASTQ file.

Ex:

```groovy
[
    id: 'FAL00870',
    strandedness: 'unstranded',
    single_end: true
]
```

\
&nbsp;

#### `reads`

Type: `path`

NextFlow input type of `path` pointing to FASTQ files on which `seqkit grep` should be run.

\
&nbsp;

#### `pattern_file`

Type: path

NextFlow input type of `path` pointing to the pattern file which has the patterns, one per line, by which FASTQ sequence ids should be searched and whose reads will be extracted.

\
&nbsp;

#### `args`

Type: Groovy String

String of optional command-line arguments to be passed to the tool. This can be mentioned in `process` scope within `withName:process_name` block using `ext.args` option within your `nextflow.config` file.

Ex:

```groovy
withName: 'SEQKIT_GREP' {
    ext.args = '--only-positive-strand'
}
```

### `output:`

___

Type: `tuple`

Outputs a tuple of metadata (`meta` from `input:`) and and filtered gzipped FASTQ file.

\
&nbsp;

#### `fastx`

Type: `path`

NextFlow output type of `path` pointing to the FASTQ format filtered gzipped file per sample (`id:`).

\
&nbsp;

#### `versions`

Type: `path`

NextFlow output type of `path` pointing to the `.yml` file storing software versions for this process.
