nextflow.enable.dsl=2

process megahit {
    tag "${id}"

    publishDir "${params.outdir}/megahit" , mode: 'copy' ,
        enabled: params.save_assembler_output ,
        pattern: "${id}/*"

    publishDir "${params.outdir}/scaffolds" , mode: 'copy' ,
        pattern: "${id}.fa"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.fa"), emit: scaffolds
    path "${id}/*", optional: true


    script:
    task_memory_GB = Math.floor(0.8 * task.memory.toGiga()) as int
    param_megahit_k = params.megahit_k == 'default' ? "" : "--k-list ${params.megahit_k}"
    input = params.single_end ? "-r \"$reads\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    param_save_assembler_output = params.save_assembler_output ? "Y" : "N"
    """
    megahit \
        $input \
        -t ${task.cpus} \
        ${param_megahit_k} \
        --memory $task_memory_GB \
        -o ${id}
    cp ${id}/final.contigs.fa ${id}.fa

    if [ "${param_save_assembler_output}" = "N" ]; then
        rm -rf ${id}
    fi
    """
}
