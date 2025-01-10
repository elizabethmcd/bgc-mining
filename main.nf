#! /usr/bin/env nextflow

// Description
// Mine biosynthetic gene clusters from bacterial genomes.

nextflow.enable.dsl=2

params.threads=16
params.outdir=null

log.info """\

MINE BGCs FROM BACTERIAL GENOMES AND ADD FUNCTIONAL ANNOTATIONS.

NOTE: YOU MUST PRE-DOWNLOAD THE ANTISMASH AND PFAM DATABASES AND PROVIDE THE PATHS
WITH --antismash_db, AND --pfam_db. THIS WORKFLOW DOES NOT SUPPORT DOWNLOADING 
DATABASES AUTOMATICALLY.
=================================================================
input_genomes                   : $params.input_genomes
genome_metadata                 : $params.genome_metadata
antismash_db                    : $params.antismash_db
pfam_db                         : $params.pfam_db
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels
// genome_fastas channel with tuple for genome name and fasta filepath
genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [file, baseName]
    }

genome_metadata = channel.fromPath(params.genome_metadata)
antismash_db_ch = channel.fromPath(params.antismash_db)
pfam_db_ch = channel.fromPath(params.pfam_db)

// workflow steps
workflow {
    // make genome STB
    all_genome_fastas_ch = genome_fastas.map{ it[0] }.collect()
    make_genome_stb(all_genome_fastas_ch)
    genome_stb_tsv = make_genome_stb.out.stb_tsv

    // predict ORFs with pyrdogial
    pyrodigal(genome_fastas)
    predicted_orfs = pyrodigal.out.predicted_orfs_gbk

    // antismash predictions
    antismash_input_ch = predicted_orfs.combine(antismash_db_ch)
    antismash(antismash_input_ch)
    antismash_gbk_files = antismash.out.gbk_results

    // bigscape on all antismash gbk_files
    all_antismash_gbk_files = antismash_gbk_files.map{ it[1] }.collect()
    run_bigscape(all_antismash_gbk_files, pfam_db_ch)
    bigscape_annotations_tsv = run_bigscape.out.bigscape_annotations_tsv

    // combine bigscape aggregate TSV with metadata
    combine_bigscape_metadata(bigscape_annotations_tsv, genome_metadata, genome_stb_tsv)

}

process make_genome_stb {
    tag "make_genome_stb"
    publishDir "${params.outdir}/genomestb", mode: 'copy'

    memory = '10 GB'
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"

    input:
    path(fasta_files)

    output:
    path("*.tsv"), emit: stb_tsv

    script:
    """
    python ${baseDir}/bin/generate_genome_stb.py ${fasta_files.join(' ')} -o genomes_stb.tsv
    """

}

process pyrodigal {
    tag "${genome_name}_pyrodigal"
    
    memory = "5 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/pyrodigal:3.4.1--py310h4b81fae_0"
    conda "envs/pyrodigal.yml"

    input:
    tuple path(fasta), val(genome_name)

    output:
    tuple val(genome_name), path("*.fna"), emit: predicted_orfs_fna
    tuple val(genome_name), path("*.faa"), emit: predicted_orfs_faa
    tuple val(genome_name), path("*.gbk"), emit: predicted_orfs_gbk

    script:
    """
    pyrodigal \\
        -i ${fasta} \\
        -f "gbk" \\
        -o "${genome_name}.gbk" \\
        -d ${genome_name}.fna \\
        -a ${genome_name}.faa
    """
}

process antismash {
    tag "${genome_name}_antismash"
    publishDir "${params.outdir}/antismash", mode: 'copy'

    memory = "20 GB"
    cpus = 4

    container "public.ecr.aws/biocontainers/antismash-lite:7.1.0--pyhdfd78af_0"
    conda "envs/antismashlite.yml"

    input:
    tuple val(genome_name), path(gbk_file), path(databases)

    output: 
    tuple val(genome_name), path("${genome_name}/*.json") , emit: json_results
    tuple val(genome_name), path("${genome_name}/*.log") , emit: log
    tuple val(genome_name), path("${genome_name}/*region*.gbk") , optional: true, emit: gbk_results

    script: 
    """
    antismash \\
        -c $task.cpus \\
        --output-dir ${genome_name} \\
        --output-basename ${genome_name} \\
        --logfile ${genome_name}/${genome_name}.log \\
        --databases $databases \\
        --genefinding-tool none \\
        ${gbk_file}
    """
}

process run_bigscape {
    tag "bigscape_all_gbks"
    publishDir "${params.outdir}/bigscape", mode: 'copy'

    memory = "36 GB"
    cpus = 10
    
    container "quay.io/biocontainers/bigscape:1.1.9--pyhdfd78af_0"
    conda "envs/bigscape.yml"

    input:
    path(gbk_files)
    path(pfam_db)

    output:
    path("*"), emit: bigscape_results
    path("bigscape_results/network_files/*/Network_Annotations_Full.tsv"), emit: bigscape_annotations_tsv

    script:
    """
    bigscape -i ./ -o bigscape_results --pfam_dir ${pfam_db} --cores ${task.cpus}
    """
}

process combine_bigscape_metadata {
    tag "combine_bigscape_metadata"
    publishDir "${params.outdir}/main_results/bgc_info", mode: 'copy'

    memory = '10 GB'
    cpus = 1

    container "public.ecr.aws/csgenetics/tidyverse:latest"
    conda "envs/tidyverse.yml"

    input:
    path(bigscape_annotations_tsv)
    path(genome_metadata)
    path(genome_stb)

    output:
    path("bgc_annotation_metadata.tsv"), emit: bgc_metadata_tsv
    path("bgc_substrate_type_counts.tsv"), emit: bgc_substrates_tsv
    path("bgc_phylo_groups_counts.tsv"), emit: bgc_phylo_groups_tsv
    
    script:
    """
    Rscript ${baseDir}/bin/combine_bgc_metadata.R ${genome_metadata} ${bigscape_annotations_tsv} ${genome_stb} bgc_annotation_metadata.tsv bgc_substrate_type_counts.tsv bgc_phylo_groups_counts.tsv
    """
}