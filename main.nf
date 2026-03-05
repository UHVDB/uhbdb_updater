#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UHVDB/uhbdb_updater
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UHVDB/uhbdb_updater
----------------------------------------------------------------------------------------
*/
// MODULES
include { SPLITBYGENUS                       } from './modules/local/splitbygenus'
include { UHBDBUPDATER as UHBDBUPDATER_LARGE } from './modules/local/uhbdbupdater'
include { UHBDBUPDATER as UHBDBUPDATER_SMALL } from './modules/local/uhbdbupdater'

// Run entry workflow
workflow {

    main:

    ch_uhbdb_dir = Channel.fromPath(params.uhbdb_dir)

    
    // MODULE: Split by genus
    SPLITBYGENUS(
        params.input,
        params.large_threshold
    )

    ch_large_input = SPLITBYGENUS.out.large
        .map { files -> files }
        .flatten()
        .map { file ->
            def genus = file.getBaseName().toString().split('_')[0]
            [ [ genus: genus ], file ]
        }
        .groupTuple(by: 0, sort: true)
        .map { meta, files -> [ meta, files[0], files[1], files[2], files[3] ] }

    UHBDBUPDATER_LARGE(
        ch_large_input,
        ch_uhbdb_dir
    )

    ch_small_input = SPLITBYGENUS.out.small
        .map { files -> files }
        .flatten()
        .map { file ->
            def genus = file.getBaseName().toString().split('_')[0]
            [ [ genus: genus ], file ]
        }
        .groupTuple(by: 0, sort: true)
        .map { meta, files -> [ meta, files[0], files[1], files[2], files[3] ] }

    UHBDBUPDATER_SMALL(
        ch_small_input,
        ch_uhbdb_dir
    )
}
