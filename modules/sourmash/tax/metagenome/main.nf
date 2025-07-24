process SOURMASH_TAX_METAGENOME {
    tag "$meta.id"
    label 'process_nano'

    module (params.enable_module ? "${params.swmodulepath}${params.fs}sourmash${params.fs}4.6.1" : null)
    conda (params.enable_conda ? "conda-forge::python bioconda::sourmash=4.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.6.1--hdfd78af_0':
        'quay.io/biocontainers/sourmash:4.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(csv), path(lineage)

    output:
    tuple val(meta), path("*.txt"), emit: txt, optional: true
    tuple val(meta), path("*.tsv"), emit: tsv, optional: true
    tuple val(meta), path("*.csv"), emit: csv, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // required defaults for the tool to run, but can be overridden
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_format = args.findAll(/(--output-format\s+[\w\,]+)\s*/).join("").replaceAll(/\,/, / --output-format /)
    args = args.replaceAll(/--output-format\s+[\w\,]+\s*/, /${output_format}/)
    """
    sourmash tax metagenome \\
        $args \\
        -g $csv \\
        --output-base $prefix \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}