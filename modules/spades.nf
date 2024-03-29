nextflow.enable.dsl=2

process metaspades {
    tag "${id}"

    publishDir "${params.outdir}/metaspades" , mode: 'copy' ,
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
    param_metaspades_k = params.metaspades_k == 'default' ? "" :  "-k ${params.metaspades_k}"
    param_metaspades_only_assembler = params.metaspades_only_assembler ? "--only-assembler" : ""
    param_save_assembler_output = params.save_assembler_output ? "Y" : "N"
    """
    spades.py \
        --meta \
        ${param_metaspades_only_assembler} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        ${param_metaspades_k} \
        --threads ${task.cpus} \
        --memory ${task_memory_GB} \
        -o ${id}
    cp ${id}/scaffolds.fasta ${id}.fa

    if [ "${param_save_assembler_output}" = "N" ]; then
        rm -rf ${id}
    fi
    """
}

process metaplasmidspades {
    tag "${id}"

    publishDir "${params.outdir}/metaplasmidspades" , mode: 'copy' ,
        enabled: params.save_assembler_output ,
        pattern: "${id}/*"
        
    publishDir "${params.outdir}/scaffolds_plasmids" , mode: 'copy' ,
        pattern: "${id}.fa"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.fa"), emit: scaffolds
    path "${id}/*", optional: true

    script:
    task_memory_GB = Math.floor(0.8 * task.memory.toGiga()) as int
    param_metaspades_k = params.metaspades_k == 'default' ? "" : "-k ${params.metaspades_k}"
    param_metaspades_only_assembler = params.metaspades_only_assembler ? "--only-assembler" : ""
    param_save_assembler_output = params.save_assembler_output ? "Y" : "N"
    """
    spades.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        --meta \
        ${param_metaspades_only_assembler} \
        --plasmid \
        ${param_metaspades_k} \
        --threads ${task.cpus} \
        --memory ${task_memory_GB} \
        -o ${id}
    cp ${id}/scaffolds.fasta ${id}.fa

    if [ "${param_save_assembler_output}" = "N" ]; then
        rm -rf ${id}
    fi
    """
}

process viralverify_db_download {

    publishDir "${params.outdir}/dbs" , mode: 'copy'

    output:
    path 'nbc_hmms.hmm', emit: viralverify_db

    script:
    """
    VIRALVERIFY_DB_URL="https://ndownloader.figshare.com/files/17904323?private_link=f897d463b31a35ad7bf0"
    curl -L \${VIRALVERIFY_DB_URL} --output nbc_hmms.hmm.gz && \
        gunzip -c nbc_hmms.hmm.gz > nbc_hmms.hmm
    rm -rf nbc_hmms.hmm.gz
    """
}

process viralverify {
    tag "${id}"

    publishDir "${params.outdir}/viralverify" , mode: 'copy' ,
        pattern: "${id}/*"

    publishDir "${params.outdir}" , mode: 'copy' ,
        saveAs: { filename ->
            if (filename == "${id}/Prediction_results_fasta/${id}_plasmid.fasta") "verified_plasmids/${id}.fa"
        }

    when:
    scaffolds.size() > 0

    input:
    tuple val(id), path(scaffolds)
    path(viralverify_db)

    output:
    path "${id}/*"
    path "${id}/Prediction_results_fasta/*"

    script:
    """
    viralverify.py \
        -f ${scaffolds} \
        -o ${id} \
        --hmm ${viralverify_db} \
        -p \
        -t ${task.cpus}
    """
}
