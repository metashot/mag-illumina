nextflow.enable.dsl=2


process deinterleave {      
    tag "${id}"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("deint*.fastq.gz"), emit: reads

    script:
    task_memory_GB = task.memory.toGiga()
    """
    reformat.sh \
        -Xmx${task_memory_GB}g \
        in=$reads \
        out1=deint_1.fastq.gz \
        out2=deint_2.fastq.gz \
        t=1
    """
}

process trim_adapters {
    tag "${id}"

    publishDir "${params.outdir}/qc/${id}" , mode: 'copy',
        pattern: "stats_adapter.txt"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("adapt*.fastq.gz"), emit: reads
    path "stats_adapter.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    output = params.single_end ? "out=adapt.fastq.gz" : "out1=adapt_1.fastq.gz out2=adapt_2.fastq.gz"
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

process remove_contaminants {
    tag "${id}"

    publishDir "${params.outdir}/qc/${id}" , mode: 'copy', 
        pattern: "stats_contaminant.txt"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("contam*.fastq.gz"), emit: reads
    path "stats_contaminant.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    output = params.single_end ? "out=contam.fastq.gz" : "out1=contam_1.fastq.gz out2=contam_2.fastq.gz"
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

process quality_filter {
    tag "${id}"

    publishDir "${params.outdir}/clean_reads" , mode: 'copy' , 
        enabled: params.save_clean

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}*.fastq.gz"), emit: reads

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    output = params.single_end ? "out=${id}.fastq.gz" : "out1=${id}_1.fastq.gz out2=${id}_2.fastq.gz"
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

process raw_reads_stats {   
    tag "${id}"

    publishDir "${params.outdir}/raw_reads_stats/${id}" , mode: 'copy'

    input:
    tuple val(id), path(reads)

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
        gcbins=auto \
        threads=${task.cpus}
    """
}

process clean_reads_stats {
    tag "${id}"

    publishDir "${params.outdir}/clean_reads_stats/${id}" , mode: 'copy'

    input:
    tuple val(id), path(reads)

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
        gcbins=auto \
        threads=${task.cpus}
    """
}

process statswrapper {      
    publishDir "${params.outdir}" , mode: 'copy'

    input:
    path scaffolds

    output:
    path 'stats_scaffolds.tsv'

    script:       
    """
    mkdir scaffolds_dir
    mv ${scaffolds} scaffolds_dir
    statswrapper.sh scaffolds_dir/* > stats_scaffolds.tsv
    """
}