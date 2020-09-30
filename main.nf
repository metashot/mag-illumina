#!/usr/bin/env nextflow

Channel
    .fromFilePairs( params.reads, size: (params.single_end || params.interleaved) ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
    .set { raw_reads_deinterleave_ch } // .into for two or more target channels

/*
 * Step 0. Deinterleave paired reads
 */
if (params.interleaved) {
    process deinterleave {      
        tag "${id}"

        input:
        tuple val(id), path(reads) from raw_reads_deinterleave_ch

        output:
        tuple val(id), path("read_*.fastq.gz") into raw_reads_stats_ch, raw_reads_adapter_ch

        script:
        task_memory_GB = task.memory.toGiga()
        
        """
        reformat.sh \
            -Xmx${task_memory_GB}g \
            in=$reads \
            out1=read_1.fastq.gz \
            out2=read_2.fastq.gz \
            t=1
        """
    }
} else {
    raw_reads_deinterleave_ch.into { raw_reads_stats_ch; raw_reads_adapter_ch }
}

/*
 * Step 1. Raw reads histograms
 */
process raw_reads_stats {   
    tag "${id}"

    publishDir "${params.outdir}/samples/${id}/raw_reads_stats" , mode: 'copy'

    input:
    tuple val(id), path(reads) from raw_reads_stats_ch

    output:
    path "*hist.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    """
    bbduk.sh \
        -Xmx${task_memory_GB}g \
        $input \
        bhist=bhist.txt \
        qhist=qhist.txt \
        gchist=gchist.txt \
        aqhist=aqhist.txt \
        lhist=lhist.txt \
        gcbins=auto
    """
}
 
/*
 * Step 2.a Remove adapters
 */
if (!params.skip_adapter) {
    process adapter {
        tag "${id}"
    
        publishDir "${params.outdir}/samples/${id}/bbduk" , mode: 'copy',
            pattern: "stats_adapter.txt"

        input:
        tuple val(id), path(reads) from raw_reads_adapter_ch

        output:
        tuple val(id), path("clean_adapter*.fastq.gz") into clean_reads_contaminant_ch
        path "stats_adapter.txt"

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean_adapter.fastq.gz" : "out1=clean_adapter_1.fastq.gz out2=clean_adapter_2.fastq.gz"
        """
        bbduk.sh \
            -Xmx${task_memory_GB}g \
            $input \
            $output \
            ref=adapters \
            ktrim=r \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe \
            tbo \
            interleaved=f \
            stats=stats_adapter.txt \
            threads=${task.cpus}
        """
    }
} else {
    raw_reads_adapter_ch.set { clean_reads_contaminant_ch }
}

/*
 * Step 2.b Remove contaminants
 */
if (!params.skip_contaminant) {
    process contaminant {
        tag "${id}"
    
        publishDir "${params.outdir}/samples/${id}/bbduk" , mode: 'copy', 
            pattern: "stats_contaminant.txt"

        input:
        tuple val(id), path(reads) from clean_reads_contaminant_ch

        output:
        tuple val(id), path("clean_contaminant*.fastq.gz") into clean_reads_quality_ch
        path "stats_contaminant.txt"

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean_contaminant.fastq.gz" : "out1=clean_contaminant_1.fastq.gz out2=clean_contaminant_2.fastq.gz"
        """
        bbduk.sh \
            -Xmx${task_memory_GB}g \
            $input \
            $output \
            ref=artifacts,phix \
            k=31 \
            hdist=1 \
            interleaved=f \
            stats=stats_contaminant.txt \
            threads=${task.cpus}
        """
    }
} else {
    clean_reads_contaminant_ch.set { clean_reads_quality_ch }
}

/*
 * Step 2.c Quality filtering/trimming and length filtering
 */
if (!params.skip_quality) {
    process quality {
        tag "${id}"

        input:
        tuple val(id), path(reads) from clean_reads_quality_ch

        output:
        tuple val(id), path("clean*.fastq.gz") into \
            clean_reads_stats_ch, clean_reads_spades_ch, clean_reads_megahit_ch, clean_reads_map_ch
  
        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean.fastq.gz" : "out1=clean_1.fastq.gz out2=clean_2.fastq.gz"
        """
        bbduk.sh \
            $input \
            $output \
            maq=10 \
            maxns=4 \
            qtrim=r \
            trimq=6 \
            mlf=0.5 \
            minlen=50 \
            interleaved=f \
            threads=${task.cpus}
            """
    }
} else {
    clean_reads_quality_ch.into {
        clean_reads_stats_ch;
        clean_reads_spades_ch;
        clean_reads_megahit_ch;
        clean_reads_map_ch 
        }
}

/*
 * Step 3. Clean reads histograms
 */

process clean_reads_stats {
    tag "${id}"

    publishDir "${params.outdir}/samples/${id}/clean_reads_stats" , mode: 'copy'

    when:
    ! (params.skip_adapter && params.skip_contaminant && params.skip_quality)

    input:
    tuple val(id), path(reads) from clean_reads_stats_ch

    output:
    path "*hist.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    """
    bbduk.sh \
        -Xmx${task_memory_GB}g \
        $input \
        bhist=bhist.txt \
        qhist=qhist.txt \
        gchist=gchist.txt \
        aqhist=aqhist.txt \
        lhist=lhist.txt \
        gcbins=auto
    """
}

/*
 * Step 4.a Assembly with Spades
 */
if (!params.single_end && !params.megahit_only) {
    process spades {
        tag "${id}"
    
        publishDir "${params.outdir}/samples/${id}" , mode: 'copy' ,
            saveAs: {filename ->
                if (filename == "spades/scaffolds.fasta") "assembly/scaffolds.fasta"
                else "$filename"
            }
    
        input:
        tuple val(id), path(reads) from clean_reads_spades_ch
    
        output:
        tuple val(id), path("spades/scaffolds.fasta") into scaffolds_spades_ch
        path "spades/spades.log"
    
        script:
        task_memory_GB = task.memory.toGiga()
        """
        spades.py \
            --meta \
            --only-assembler \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            --threads ${task.cpus} \
            --memory ${task_memory_GB} \
            -o spades
        """
    }
} else {
    scaffolds_spades_ch = Channel.empty()
}

/*
 * Step 4.b Assembly with Megahit
 */

if (params.single_end || params.megahit_only) {
    process megahit {
        tag "${id}"

        publishDir "${params.outdir}/samples/${id}" , mode: 'copy' ,
            saveAs: {filename -> if (filename == "megahit/final.contigs.fa") "assembly/scaffolds.fasta"}

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

    publishDir "${params.outdir}/samples/${id}/assembly" , mode: 'copy'

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
    per_thread_memory_GB = Math.floor(task.memory.toGiga() / task.cpus) as int
    per_thread_memory_GB = per_thread_memory_GB > 1 ? per_thread_memory_GB : 1
    """
    samtools sort -@ ${task.cpus} -m ${per_thread_memory_GB}G -o map.bam $sam
    """
}

scaffolds_metabat2_ch.join(map_metabat2_ch).set{ merged_metabat2_ch }

/*
 * Step 6.c Binning with Metabat2
 */
process metabat2 {
    tag "${id}"

    publishDir "${params.outdir}/samples/${id}" , mode: 'copy' ,
        saveAs: {filename ->
            if (filename == "metabat2_depth.txt") "assembly/depth.txt"
            else if (filename == "metabat2.log") "metabat2/log.txt"
            else "$filename"
        }

    input:
    tuple val(id), path(scaffolds), path(bam) from merged_metabat2_ch
    
    output:
    path "bins/*.bin*" optional true
    path "metabat2_depth.txt"
    path "metabat2.log"
    
    script:
    """
    mkdir -p bins
    jgi_summarize_bam_contig_depths --outputDepth metabat2_depth.txt $bam
    metabat2 \
        -i $scaffolds \
        -a metabat2_depth.txt \
        -o bins/${id}.bin \
        -v \
        -t ${task.cpus} \
        -m 2500 &> metabat2.log
    """
}
