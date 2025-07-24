// Help text for gen_sim_abn_table.py (gsat) within CPIPES.

def gsatpyHelp(params) {

    Map tool = [:]
    Map toolspecs = [:]
    tool.text = [:]
    tool.helpparams = [:]

    toolspecs = [
        'gsatpy_run': [
            clihelp: 'Run the gen_sim_abn_table.py script. Default: ' +
                (params.gsatpy_run ?: false),
            cliflag: null,
            clivalue: null
        ],
        'gsatpy_header': [
            clihelp: 'Does the taxonomic summary result files have ' +
                'a header line. ' +
                " Default: ${params.gsatpy_header}",
            cliflag: '-header',
            clivalue: (params.gsatpy_header ? ' ' : '')
        ]
    ]

    toolspecs.each {
        k, v -> tool.text['--' + k] = "${v.clihelp}"
        tool.helpparams[k] = [ cliflag: "${v.cliflag}", clivalue: v.clivalue ]
    }

    return tool
}