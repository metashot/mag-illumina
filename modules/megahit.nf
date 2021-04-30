nextflow.enable.dsl=2

process megahit {
    tag "${id}"

    publishDir "${params.outdir}/megahit" , mode: 'copy' ,
        pattern: "${id}/*"

    publishDir "${params.outdir}/scaffolds" , mode: 'copy' ,
        pattern: "${id}.fa"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.fa"), emit: scaffolds
    path "${id}/*"


    script:
    task_memory_GB = task.memory.toGiga()
    param_megahit_k = params.megahit_k == 'default' ? "" :  "-k-list ${params.megahit_k}"
    input = params.single_end ? "-r \"$reads\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    megahit \
        $input \
        -t ${task.cpus} \
        ${param_megahit_k} \
        --memory $task_memory_GB \
        -o ${id}
    cp ${id}/final.contigs.fa ${id}.fa
    """
}
