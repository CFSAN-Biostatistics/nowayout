// Help text for `gen_salmon_tph_and_krona_tsv.py` (gsalkronapy) within CPIPES.

def gsalkronapyHelp(params) {

    Map tool = [:]
    Map toolspecs = [:]
    tool.text = [:]
    tool.helpparams = [:]

    toolspecs = [
        'gsalkronapy_run': [
            clihelp: 'Run the `gen_salmon_tph_and_krona_tsv.py` script. Default: ' +
                (params.gsalkronapy_run ?: false),
            cliflag: null,
            clivalue: null
        ],
        'gsalkronapy_sf': [
            clihelp: 'Set the scaling factor by which TPM values ' +
                'are scaled down.' +
                " Default: ${params.gsalkronapy_sf}",
            cliflag: '-sf',
            clivalue: (params.gsalkronapy_sf ?: '')
        ],
        'gsalkronapy_smres_suffix': [
            clihelp: 'Find the `sourmash gather` result files ' +
                'ending in this suffix.' +
                " Default: ${params.gsalkronapy_smres_suffix}",
            cliflag: '-smres-suffix',
            clivalue: (params.gsalkronapy_smres_suffix ?: '')
        ],
        'gsalkronapy_failed_suffix': [
            clihelp: 'Find the sample names which failed classification stored ' +
                'inside the files ending in this suffix.' +
                " Default: ${params.gsalkronapy_failed_suffix}",
            cliflag: '-failed-suffix',
            clivalue: (params.gsalkronapy_failed_suffix ?: '')
        ],
        'gsalkronapy_num_lin_cols': [
            clihelp: 'Number of columns expected in the lineages CSV file. ' +
                " Default: ${params.gsalkronapy_num_lin_cols}",
            cliflag: '-num-lin-cols',
            clivalue: (params.gsalkronapy_num_lin_cols ?: '')
        ],
        'gsalkronapy_lin_regex': [
            clihelp: 'Number of columns expected in the lineages CSV file. ' +
                " Default: ${params.gsalkronapy_num_lin_cols}",
            cliflag: '-num-lin-cols',
            clivalue: (params.gsalkronapy_num_lin_cols ?: '')
        ]
    ]

    toolspecs.each {
        k, v -> tool.text['--' + k] = "${v.clihelp}"
        tool.helpparams[k] = [ cliflag: "${v.cliflag}", clivalue: v.clivalue ]
    }

    return tool
}