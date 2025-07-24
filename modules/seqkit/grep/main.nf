process SEQKIT_GREP {
    tag "$meta.id"
    label 'process_low'

    module (params.enable_module ? "${params.swmodulepath}${params.fs}seqkit${params.fs}2.2.0" : null)
    conda (params.enable_conda ? "bioconda::seqkit=2.2.0 conda-forge::sed=4.7 conda-forge::coreutils" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0':
        'quay.io/biocontainers/seqkit:2.1.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads), path(pattern_file)

    output:
    tuple val(meta), path("*.gz"), emit: fastx
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num_read_files = reads.toList().size()
    def extension = "fastq"
    if ("$reads" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz/) {
        extension = "fasta"
    }

    if (meta.single_end || num_read_files == 1) {
        """
        pattern_file_contents=\$(sed '1!d' $pattern_file)
        if [ "\$pattern_file_contents" != "DuMmY" ]; then
            cut -f1 -d " " $pattern_file > ${prefix}.seqids.txt
            additional_args="-f ${prefix}.seqids.txt $args"
        else
            additional_args="$args"
        fi

        seqkit \\
            grep \\
            -j $task.cpus \\
            -o ${prefix}.seqkit-grep.${extension}.gz \\
            \$additional_args \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
    } else {
        """
        pattern_file_contents=\$(sed '1!d' $pattern_file)
        if [ "\$pattern_file_contents" != "DuMmY" ]; then
            additional_args="-f $pattern_file $args"
        else
            additional_args="$args"
        fi

        seqkit \\
            grep \\
            -j $task.cpus \\
            -o ${prefix}.R1.seqkit-grep.${extension}.gz \\
            \$additional_args \\
            ${reads[0]}
        
        seqkit \\
            grep \\
            -j $task.cpus \\
            -o ${prefix}.R2.seqkit-grep.${extension}.gz \\
            \$additional_args \\
            ${reads[1]}

        seqkit \\
            pair \\
            -j $task.cpus \\
            -1 ${prefix}.R1.seqkit-grep.${extension}.gz \\
            -2 ${prefix}.R2.seqkit-grep.${extension}.gz

        rm ${prefix}.R1.seqkit-grep.${extension}.gz
        rm ${prefix}.R2.seqkit-grep.${extension}.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
    }
}
