// Define any required imports for this specific workflow
import java.nio.file.Paths
import java.util.zip.GZIPInputStream
import java.io.FileInputStream
import nextflow.file.FileHelper


// Include any necessary methods
include { \
    summaryOfParams; stopNow; fastqEntryPointHelp; sendMail; \
    addPadding; wrapUpHelp           } from "${params.routines}"
include { fastpHelp                  } from "${params.toolshelp}${params.fs}fastp"
include { kmaalignHelp               } from "${params.toolshelp}${params.fs}kmaalign"
include { seqkitgrepHelp             } from "${params.toolshelp}${params.fs}seqkitgrep"
include { salmonidxHelp              } from "${params.toolshelp}${params.fs}salmonidx"
include { sourmashsketchHelp         } from "${params.toolshelp}${params.fs}sourmashsketch"
include { sourmashgatherHelp         } from "${params.toolshelp}${params.fs}sourmashgather"
include { sfhpyHelp                  } from "${params.toolshelp}${params.fs}sfhpy"
include { gsalkronapyHelp            } from "${params.toolshelp}${params.fs}gsalkronapy"
include { kronaktimporttextHelp      } from "${params.toolshelp}${params.fs}kronaktimporttext"

// Exit if help requested before any subworkflows
if (params.help) {
    log.info help()
    exit 0
}


// Include any necessary modules and subworkflows
include { PROCESS_FASTQ           } from "${params.subworkflows}${params.fs}process_fastq"
include { FASTP                   } from "${params.modules}${params.fs}fastp${params.fs}main"
include { KMA_ALIGN               } from "${params.modules}${params.fs}kma${params.fs}align${params.fs}main"
include { OTF_GENOME              } from "${params.modules}${params.fs}otf_genome${params.fs}main"
include { SEQKIT_GREP             } from "${params.modules}${params.fs}seqkit${params.fs}grep${params.fs}main"
include { SALMON_INDEX            } from "${params.modules}${params.fs}salmon${params.fs}index${params.fs}main"
include { SALMON_QUANT            } from "${params.modules}${params.fs}salmon${params.fs}quant${params.fs}main"
include { SOURMASH_SKETCH         } from "${params.modules}${params.fs}sourmash${params.fs}sketch${params.fs}main"
include { SOURMASH_SKETCH \
    as REDUCE_DB_IDX              } from "${params.modules}${params.fs}sourmash${params.fs}sketch${params.fs}main"
include { SOURMASH_GATHER         } from "${params.modules}${params.fs}sourmash${params.fs}gather${params.fs}main"
include { NOWAYOUT_RESULTS        } from "${params.modules}${params.fs}nowayout_results${params.fs}main"
include { KRONA_KTIMPORTTEXT      } from "${params.modules}${params.fs}krona${params.fs}ktimporttext${params.fs}main"
include { DUMP_SOFTWARE_VERSIONS  } from "${params.modules}${params.fs}custom${params.fs}dump_software_versions${params.fs}main"
include { MULTIQC                 } from "${params.modules}${params.fs}multiqc${params.fs}main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUTS AND ANY CHECKS FOR THE BETTERCALLSAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def reads_platform = 0
reads_platform += (params.input ? 1 : 0)

if (reads_platform < 1 || reads_platform == 0) {
    stopNow("Please mention at least one absolute path to input folder which contains\n" +
            "FASTQ files sequenced using the --input option.\n" +
        "Ex: --input (Illumina or Generic short reads in FASTQ format)")
}

params.fastp_adapter_fasta ? checkMetadataExists(params.fastp_adapter_fasta, 'Adapter sequences FASTA') : null
checkMetadataExists(params.lineages_csv, 'Lineages CSV')
checkMetadataExists(params.kmaalign_idx, 'KMA Indices')
checkMetadataExists(params.ref_fna, 'FASTA reference')

ch_sourmash_lin = file( params.lineages_csv )


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN THE BETTERCALLSAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NOWAYOUT {
    main:
        log.info summaryOfParams()

        PROCESS_FASTQ()

        PROCESS_FASTQ.out.versions
            .set { software_versions }

        PROCESS_FASTQ.out.processed_reads
            .set { ch_processed_reads }

        ch_processed_reads
            .map { meta, fastq ->
                meta.get_kma_hit_accs = true
                meta.salmon_decoys = params.dummyfile
                meta.salmon_lib_type = (params.salmonalign_libtype ?: false)
                meta.kma_t_db = params.kmaalign_idx
                [ meta, fastq ]
            }
            .filter { meta, fastq ->
                fq_file = ( fastq.getClass().toString() =~ /ArrayList/ ? fastq : [ fastq ] )
                fq_gzip = new GZIPInputStream( new FileInputStream( fq_file[0].toAbsolutePath().toString() ) )
                fq_gzip.read() != -1
            }
            .set { ch_processed_reads }

        FASTP( ch_processed_reads )

        FASTP.out.json
            .map { meta, json ->
                json
            }
            .collect()
            .set { ch_multiqc }

        KMA_ALIGN(
            FASTP.out.passed_reads
                .map { meta, fastq ->
                    [meta, fastq, []]
                }
        )

        OTF_GENOME(
            KMA_ALIGN.out.hits
                .join(KMA_ALIGN.out.frags)
        )

        OTF_GENOME.out.reads_extracted
            .filter { meta, fasta ->
                fa_file = ( fasta.getClass().toString() =~ /ArrayList/ ? fasta : [ fasta ] )
                fa_gzip = new GZIPInputStream( new FileInputStream( fa_file[0].toAbsolutePath().toString() ) )
                fa_gzip.read() != -1
            }
            .set { ch_mito_aln_reads }

        SEQKIT_GREP(
            KMA_ALIGN.out.hits
                .filter { meta, mapped_refs ->
                    patterns = file( mapped_refs )
                    patterns.size() > 0
                }
                .map { meta, mapped_refs ->
                    [meta, params.ref_fna, mapped_refs]
                }
        )

        SALMON_INDEX( SEQKIT_GREP.out.fastx )

        SALMON_QUANT(
            ch_mito_aln_reads
                .join( SALMON_INDEX.out.idx )
        )

        REDUCE_DB_IDX(
            SEQKIT_GREP.out.fastx,
            true,
            false,
            'db'
        )

        SOURMASH_SKETCH(
            ch_mito_aln_reads,
            false,
            false,
            'query'
        )

        SOURMASH_GATHER(
            SOURMASH_SKETCH.out.signatures
                .join( REDUCE_DB_IDX.out.signatures ),
                [], [], [], []
        )

        // SOURMASH_TAX_METAGENOME(
        //     SOURMASH_GATHER.out.result
        //         .groupTuple(by: [0])
        //         .map { meta, csv ->
        //             [ meta, csv, ch_sourmash_lin ]
        //         }
        // )

        // SOURMASH_TAX_METAGENOME.out.csv
        //     .map { meta, csv ->
        //         csv
        //     }
        //     .set { ch_lin_csv }

        // SOURMASH_TAX_METAGENOME.out.tsv
        //     .tap { ch_lin_krona }
        //     .map { meta, tsv ->
        //         tsv
        //     }
        //     .tap { ch_lin_tsv }

        SOURMASH_GATHER.out.result
            .groupTuple(by: [0])
            .map { meta, csv ->
                [ csv ]
            }
            .concat(
                SALMON_QUANT.out.results
                    .map { meta, salmon_res ->
                        [ salmon_res ]
                    }
            )
            .concat(
                SOURMASH_GATHER.out.failed
                    .map { meta, failed ->
                        [ failed ]
                    }
            )
            .concat( OTF_GENOME.out.failed )
            .collect()
            .flatten()
            .collect()
            .set { ch_gene_abn }
        
        NOWAYOUT_RESULTS( ch_gene_abn, ch_sourmash_lin )

        NOWAYOUT_RESULTS.out.tsv
            .flatten()
            .filter { tsv -> tsv.toString() =~ /.*${params.krona_res_suffix}$/ }
            .map { tsv ->
                    meta = [:]
                    meta.id = "${params.cfsanpipename}_${params.pipeline}_krona"
                    [ meta, tsv ]
            }
            .groupTuple(by: [0])
            .set { ch_lin_krona }

        // ch_lin_tsv
        //     .mix( ch_lin_csv )
        //     .collect()
        //     .set { ch_lin_summary }

        // SOURMASH_TAX_METAGENOME.out.txt
        //     .map { meta, txt ->
        //         txt
        //     }
        //     .collect()
        //     .set { ch_lin_kreport }

        // NOWAYOUT_RESULTS(
        //     ch_lin_summary
        //         .concat( SOURMASH_GATHER.out.failed )
        //         .concat( OTF_GENOME.out.failed )
        //         .collect()
        // )

        KRONA_KTIMPORTTEXT( ch_lin_krona )
        
        DUMP_SOFTWARE_VERSIONS(
            software_versions
                .mix (
                    FASTP.out.versions,
                    KMA_ALIGN.out.versions,
                    SEQKIT_GREP.out.versions,
                    REDUCE_DB_IDX.out.versions,
                    SOURMASH_SKETCH.out.versions,
                    SOURMASH_GATHER.out.versions,
                    SALMON_INDEX.out.versions,
                    SALMON_QUANT.out.versions,
                    NOWAYOUT_RESULTS.out.versions,
                    KRONA_KTIMPORTTEXT.out.versions
                )
                .unique()
                .collectFile(name: 'collected_versions.yml')
        )

        DUMP_SOFTWARE_VERSIONS.out.mqc_yml
            .concat(
                ch_multiqc,
                NOWAYOUT_RESULTS.out.mqc_yml
            )
            .collect()
            .flatten()
            .collect()
            .set { ch_multiqc }

        MULTIQC( ch_multiqc )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ON COMPLETE, SHOW GORY DETAILS OF ALL PARAMS WHICH WILL BE HELPFUL TO DEBUG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success) {
        sendMail()
    }
}

workflow.onError {
    sendMail()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHOD TO CHECK METADATA EXISTENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkMetadataExists(file_path, msg) {
    file_path_obj = file( file_path )
    
    if (msg.toString().find(/(?i)KMA/)) {
        if (!file_path_obj.parent.exists()) {
            stopNow("Please check if your ${msg}\n" +
                "[ ${file_path} ]\nexists and that the files are not of size 0.")
        }
        
        // Check if db files within parent path are empty.
        file_path_obj.parent.eachFileRecurse {
            if (it.size() == 0) {
                stopNow("For ${msg}, within\n" +
                "[ ${file_path} ],\nthe following file is of size 0: ${it.name}")
            }
        }

    }
    else if (!file_path_obj.exists() || file_path_obj.size() == 0) {
        stopNow("Please check if your ${msg} file\n" +
            "[ ${file_path} ]\nexists and is not of size 0.")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP TEXT METHODS FOR BETTERCALLSAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def help() {

    Map helptext = [:]

    helptext.putAll (
        fastqEntryPointHelp() +
        fastpHelp(params).text +
        kmaalignHelp(params).text +
        seqkitgrepHelp(params).text +
        salmonidxHelp(params).text +
        sourmashsketchHelp(params).text +
        sourmashgatherHelp(params).text +
        sfhpyHelp(params).text +
        gsalkronapyHelp(params).text +
        kronaktimporttextHelp(params).text +
        wrapUpHelp()
    )

    return addPadding(helptext)
}
