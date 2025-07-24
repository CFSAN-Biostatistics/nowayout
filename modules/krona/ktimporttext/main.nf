process KRONA_KTIMPORTTEXT {
    tag "$meta.id"
    label 'process_nano'

    module (params.enable_module ? "${params.swmodulepath}${params.fs}krona${params.fs}2.8.1" : null)
    conda (params.enable_conda ? "conda-forge::curl bioconda::krona=2.8.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1':
        'quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1' }"

    input:
        tuple val(meta), path(report)

    output:
        tuple val(meta), path ('*.html'), emit: html
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def krona_suffix = params.krona_res_suffix ?: '.krona.tsv'
        def reports = report.collect {
            it = it.toString() + ',' + it.toString().replaceAll(/(.*)${krona_suffix}$/, /$1/)
        }.sort().join(' ')
        """
        ktImportText  \\
            $args \\
            -o ${prefix}.html \\
            $reports

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
        END_VERSIONS
        """
}