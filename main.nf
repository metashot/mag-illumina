#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { deinterleave; trim_adapters; remove_contaminants; quality_filter; raw_reads_stats; clean_reads_stats; statswrapper } from './modules/bbtools'
include { metaspades; metaplasmidspades; viralverify; viralverify_db_download} from './modules/spades'
include { megahit } from './modules/megahit'
include { map; sam2bam; metabat2 } from './modules/metabat'


workflow {
    Channel
        .fromFilePairs( params.reads, size: (params.single_end || params.interleaved) ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
        .set { reads_ch }

    if (params.interleaved) {
        deinterleave(reads_ch)
        deint_ch = deinterleave.out.reads
    } else {
        deint_ch = reads_ch
    }

    raw_reads_stats(deint_ch)

    if (!params.skip_cleaning) {
        trim_adapters(deint_ch)
        remove_contaminants(trim_adapters.out.reads)
        quality_filter(remove_contaminants.out.reads)
        clean_ch = quality_filter.out.reads
        clean_reads_stats(clean_ch)
    } else {
        clean_ch = reads_ch
    }

    if (params.run_metaplasmidspades && !params.single_end) {
        metaplasmidspades(clean_ch)

        if (params.viralverify_db == 'none') {
            viralverify_db_download()
            viralverify_db = viralverify_db_download.out.viralverify_db
        }
        else {
            viralverify_db = file(params.viralverify_db, checkIfExists: true)
        }

        viralverify(metaplasmidspades.out.scaffolds, viralverify_db)
    }

    if (!params.skip_assembly){
        if ((params.assembler == "metaspades") && !params.single_end) {
            metaspades(clean_ch)
            scaffolds_ch = metaspades.out.scaffolds
        }
        else {
            megahit(clean_ch)
            scaffolds_ch = megahit.out.scaffolds
        }

        statswrapper(scaffolds_ch.map { row -> row[1] }.collect())

        map(clean_ch.join(scaffolds_ch))
        sam2bam(map.out.sam)
        metabat2(scaffolds_ch.join(sam2bam.out.bam))
    }   
}
