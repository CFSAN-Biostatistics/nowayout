process SOURMASH_SKETCH {
    tag "$meta.id"
    label 'process_nano'

    module (params.enable_module ? "${params.swmodulepath}${params.fs}sourmash${params.fs}4.6.1" : null)
    conda (params.enable_conda ? "conda-forge::python bioconda::sourmash=4.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.6.1--hdfd78af_0':
        'quay.io/biocontainers/sourmash:4.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(sequence)
    val singleton
    val merge
    val db_or_query

    output:
    tuple val(meta), path("*.{query,db}.sig"), emit: signatures
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // required defaults for the tool to run, but can be overridden
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_sig = merge ? "--merge ${meta.id}" : ''
    def singleton = singleton ? '--singleton' : ''
    """
    sourmash sketch \\
        $args \\
        $merge_sig \\
        $singleton \\
        --output "${prefix}.${db_or_query}.sig" \\
        $sequence

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sourmash: \$(echo \$(sourmash --version 2>&1) | sed 's/^sourmash //' )
    END_VERSIONS
    """
}