nextflow.enable.dsl=2

process map {
    tag "${id}"

    input:
    tuple val(id), path(reads), path(scaffolds)
    
    output:
    tuple val(id), path("map.sam"), emit: sam
    
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

process sam2bam {
    tag "${id}"

    input:
    tuple val(id), path(sam)
    
    output:
    tuple val(id), path("map.bam"), emit: bam
    
    script:
    thread_mem_GB = Math.floor(task.memory.toGiga() / task.cpus) - 2) as int
    thread_mem_GB_1 = thread_mem_GB > 1 ? thread_mem_GB : 1
    """
    samtools sort -@ ${task.cpus} -m ${thread_mem_GB_1}G -o map.bam $sam
    """
}

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
    tuple val(id), path(scaffolds), path(bam)
    
    output:
    path "bins/*.fa", optional: true
    path "unbinned/*.fa", optional: true
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
