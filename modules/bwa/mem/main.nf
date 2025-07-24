process BWA_MEM {
    tag "$meta.id"
    label 'process_micro'

    module (params.enable_module ? "${params.swmodulepath}${params.fs}bwa${params.fs}0.7.17" : null)
    conda (params.enable_conda ? "bioconda::bwa=0.7.17 conda-forge::perl" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--he4a0461_11' :
        'quay.io/biocontainers/bwa:0.7.17--he4a0461_11' }"

    input:
        tuple val(meta), path(reads), path(index)
        val index2

    output:
        tuple val(meta), path("*.sam"), emit: aligned_sam
        path  "versions.yml"          , emit: versions

    when:
        

    script:
        def args   = task.ext.args ?: ''
        def args2  = task.ext.args2 ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def this_index = (index ?: index2)
        """

        if [ "${params.fq_single_end}" = "false" ]; then
            bwa mem \\
                $args \\
                -t $task.cpus \\
                $this_index \\
                ${reads[0]} ${reads[1]} > ${prefix}.aligned.sam
        else
            bwa mem \\
                $args \\
                -t $task.cpus \\
                -a \\
                $this_index \\
                $reads > ${prefix}.aligned.sam

        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        END_VERSIONS
        """
}