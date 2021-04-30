#!/usr/bin/env nextflow

extflow.enable.dsl=2

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
    }

    if (params.run_metaplasmidspades) {
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
}






/*
 * Step 4.b Assembly with Megahit
 */

if (params.single_end || params.megahit_only) {
    process megahit {
        tag "${id}"

        publishDir "${params.outdir}" , mode: 'copy' ,
            saveAs: { filename -> 
                if (filename == "megahit/final.contigs.fa") "scaffolds/${id}.scaffolds.fa"
            }

        input:
        tuple val(id), path(reads) from clean_reads_megahit_ch

        output:
        tuple val(id), path("megahit/final.contigs.fa") into scaffolds_megahit_ch

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "-r \"$reads\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
        """
        megahit \
            $input \
            -t ${task.cpus} \
            --k-min 27 \
            --k-max 99 \
            --k-step 14 \
            --kmin-1pass \
            --memory $task_memory_GB \
            -o megahit
        """
    }
} else {
    scaffolds_megahit_ch = Channel.empty()
}

scaffolds_spades_ch
    .mix(scaffolds_megahit_ch)
    .into { scaffolds_stats_ch; scaffolds_map_ch; scaffolds_metabat2_ch }

/*
 * Step 5. Scaffold statistics
 */
process assembly_stats {
    tag "${id}"

    publishDir "${params.outdir}/assembly_stats/${id}" , mode: 'copy'

    input:
    tuple val(id), path(scaffolds) from scaffolds_stats_ch

    output:
    path "stats.txt"

    script:
    """
    stats.sh \
        in=$scaffolds \
        out=stats.txt
    """
}

clean_reads_map_ch.join(scaffolds_map_ch).set{ merged_map_ch } 

/*
 * Step 6.a Map reads to contigs
 */
process map {
    tag "${id}"

    input:
    tuple val(id), path(reads), path(scaffolds) from merged_map_ch
    
    output:
    tuple val(id), path("map.sam") into map_sam2bam_ch
    
    script:
    input = params.single_end ? "-U \"$reads\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    mkdir -p bowtie_db

    bowtie2-build \
        --threads ${task.cpus} \
        $scaffolds \
        bowtie_db/scaffolds

    bowtie2 \
        -x bowtie_db/scaffolds \
        $input \
        --threads ${task.cpus} \
        -S map.sam
    """
}

/*
 * Step 6.b SAM2BAM
 */
process sam2bam {
    tag "${id}"

    input:
    tuple val(id), path(sam) from map_sam2bam_ch
    
    output:
    tuple val(id), path("map.bam") into map_metabat2_ch
    
    script:
    thread_mem_GB = Math.floor(task.memory.toGiga() / task.cpus) as int
    thread_mem_GB_1 = thread_mem_GB > 1 ? thread_mem_GB : 1
    """
    samtools sort -@ ${task.cpus} -m ${thread_mem_GB_1}G -o map.bam $sam
    """
}

scaffolds_metabat2_ch.join(map_metabat2_ch).set{ merged_metabat2_ch }

/*
 * Step 6.c Binning with Metabat2
 */
process metabat2 {
    tag "${id}"

    publishDir "${params.outdir}" , mode: 'copy' ,
        pattern: "{bins, unbinned}/*.fa"

    publishDir "${params.outdir}" , mode: 'copy' ,
        saveAs: {filename ->
            if (filename == "metabat2_depth.txt") "metabat2/${id}/depth.txt"
            else if (filename == "metabat2.log") "metabat2/${id}/log.txt"
        }

    input:
    tuple val(id), path(scaffolds), path(bam) from merged_metabat2_ch
    
    output:
    path "bins/*.fa"  optional true
    path "unbinned/*.fa"  optional true
    path "metabat2_depth.txt"
    path "metabat2.log"
    
    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth metabat2_depth.txt $bam

    mkdir -p bins
    metabat2 \
        -i $scaffolds \
        -a metabat2_depth.txt \
        -o bins/${id}.bin \
        -v \
        -t ${task.cpus} \
        --seed 1 \
        --unbinned \
        -m ${params.min_contig_size} &> metabat2.log

    mkdir -p unbinned
    if [ -f bins/${id}.bin.unbinned.fa ]; then
        mv bins/${id}.bin.unbinned.fa unbinned/${id}.unbinned.fa
    fi
    """
}
