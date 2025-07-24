process NOWAYOUT_RESULTS {
    tag "nowayout aggregate"
    label "process_pico"

    module (params.enable_module ? "${params.swmodulepath}${params.fs}python${params.fs}3.8.1" : null)
    conda (params.enable_conda ? 'conda-forge::python=3.11 conda-forge::spectra conda-forge::lzstring conda-forge::imp bioconda::multiqc=1.19' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
        path pass_and_fail_rel_abn_files
        path lineage_csv

    output:
        path '*.tblsum.txt', emit: mqc_txt, optional: true
        path '*_mqc.json'  , emit: mqc_json, optional: true
        path '*_mqc.yml'   , emit: mqc_yml, optional: true
        path '*.tsv'       , emit: tsv, optional: true
        path 'versions.yml', emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        gen_salmon_tph_and_krona_tsv.py \\
            $args \\
            -sal "." \\
            -smres "." \\
            -lin $lineage_csv

        create_mqc_data_table.py \\
            "nowayout" "The results shown here are <code>salmon quant</code> TPM values scaled down by a factor of ${params.gsalkronapy_sf}."

        create_mqc_data_table.py \\
            "nowayout_indiv_reads_mapped" "The results shown here are the number of reads mapped (post threshold filters) per taxon to the <code>nowayout</code>'s custom <code>${params.db_mode}</code> database for each sample."

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$( python --version | sed 's/Python //g' )
        END_VERSIONS
        """
}