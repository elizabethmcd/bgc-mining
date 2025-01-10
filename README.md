# Mining biosynthetic gene clusters from bacterial genomes

This repository contains the `bgc-mining` Nextflow workflow for mining biosynthetic gene clusters from bacterial genomes.

The main input is a directory of bacterial genomes in FASTA format. The workflow will mine for biosynthetic gene clusters (BGCs) using antiSMASH and summarize the BGC clusters with Bigscape. Note that this worfklow is configured for custom purposes such as inputting a specific genome metadata TSV file for joining metadata with main results files. You can see an example of the required input metadata file in the `test_data/metadata` directory.

## Workflow Usage
The pipeline can be run with either conda or docker using the `-profile` flag. All input genomes in the input directory should end in `.fa`. 

Importantly, the workflow does not handle automatic downloading of the antiSMASH and PFAM databases. You must pre-download these databases and provide the paths with `--antismash_db` and `--pfam_db`.
```
nextflow run main.nf \\
    --input_genomes <INPUT_DIRECTORY> \\
    --outdir <OUTPUT_DIRECTORY> \\
    --antismash_db <ANTISMASH_DB_PATH> \\
    --pfam_db <PFAM_DB_PATH> \\
    --genome_metadata <GENOME_METADATA_PATH> \\
    -profile <docker|conda>
```

## Databases
This workflow requires pre-downloading the antiSMASH and PFAM databases. The PFAM database requried for Bigscape is already available within the antiSMASH download. To download the antiSMASH database, it's recommended that you use the database downloader script available in the command-line version of antiSMASH. See [here](https://docs.antismash.secondarymetabolites.org/install/) for documentation for how to install antiSMASH and use the download helper script. 