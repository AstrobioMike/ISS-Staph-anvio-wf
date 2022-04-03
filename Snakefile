####################################################################
## Snakefile for Lee and O'Rourke et al. ISS Staph anvio workflow ##
####################################################################

import os

configfile: "config.yaml"


########################################
############# General Info #############
########################################

"""
See the corresponding 'config.yaml' file for general use information.
Variables that may need to be adjusted should be changed there.
"""

########################################
#### Reading samples file into list ####
########################################

genome_IDs_file = config["genome_IDs_file"]
genome_ID_list = [line.strip() for line in open(genome_IDs_file)]



########################################
######## Setting up directories ########
########################################


dirs_to_create = [
                  config["input_anvio_files_dir"], 
                  config["bam_files_dir"], 
                  config["bowtie2_indexes_dir"], 
                  config["contigs_dbs_dir"], 
                  config["profile_dbs_dir"]
                  ]

for dir in dirs_to_create:
    try:
        os.mkdir(dir)
    except:
        pass


########################################
############# Rules start ##############
########################################


rule all:
    input:
        anvio_input_files = expand(os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa"), genome = genome_ID_list),
        bowtie2_indexes = expand(os.path.join(config["bowtie2_indexes_dir"], "{genome}.1.bt2"), genome = genome_ID_list),
        bam_files = expand(os.path.join(config["bam_files_dir"], "{genome}.bam"), genome = genome_ID_list),
        contigs_dbs = expand(os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db"), genome = genome_ID_list),
        profile_dbs = expand(os.path.join(config["profile_dbs_dir"], "{genome}-profile"), genome = genome_ID_list),
        external_genomes_file = "external-genomes.tsv",
        genomes_db = "our-isolates-GENOMES.db",
        pan = config["anvio_pan_output_dir"],
        

rule genbank_to_anvio:
    input:
        genbank = os.path.join(config["genbank_files_dir"], "{genome}.gb")
    output:
        fasta = os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa"),
        gene_calls = os.path.join(config["input_anvio_files_dir"], "{genome}-external-gene-calls.txt"),
        gene_functions = os.path.join(config["input_anvio_files_dir"], "{genome}-external-functions.txt")
    params:
        output_prefix = os.path.join(config["input_anvio_files_dir"], "{genome}")
    log:
        config["logs_dir"] + "genbank_to_anvio-{genome}.log"
    shell:
        """
        anvi-script-process-genbank -i {input.genbank} -O {params.output_prefix} > {log} 2>&1
        """


rule make_bowtie2_index:
    input:
        fasta = os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa")
    output:
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.1.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.2.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.3.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.4.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.rev.1.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.rev.2.bt2")
    params:
        output_prefix = os.path.join(config["bowtie2_indexes_dir"], "{genome}"),
        num_threads = config["mapping_threads"]
    log:
        config["logs_dir"] + "make_bowtie2_index-{genome}.log"
    shell:
        """
        bowtie2-build --threads {params.num_threads} {input.fasta} {params.output_prefix} > {log} 2>&1
        """


rule mapping:
    input:
        index_built_trigger = os.path.join(config["bowtie2_indexes_dir"], "{genome}.1.bt2"),
        R1 = os.path.join(config["reads_dir"] + "{genome}-err-corr-R1.fq.gz"),
        R2 = os.path.join(config["reads_dir"] + "{genome}-err-corr-R2.fq.gz")
    output:
        bam = os.path.join(config["bam_files_dir"] + "{genome}.bam")
    params:
        index = os.path.join(config["bowtie2_indexes_dir"], "{genome}"),
        num_threads = config["mapping_threads"]
    log:
        os.path.join(config["bam_files_dir"] + "{genome}-mapping-info.txt")
    shell:
        """
        bowtie2 -q --threads {params.num_threads} -x {params.index} -1 {input.R1} -2 {input.R2} --no-unal 2> {log} | samtools view -b | samtools sort -@ {params.num_threads} > {output.bam}
        samtools index -@ {params.num_threads} {output.bam}
        """


rule make_and_annotate_contigs_db:
    input:
        fasta = os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa"),
        gene_calls = os.path.join(config["input_anvio_files_dir"], "{genome}-external-gene-calls.txt"),
        gene_functions = os.path.join(config["input_anvio_files_dir"], "{genome}-external-functions.txt"),
    output:
        contigs_db = os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db")
    params:
        name = "{genome}",
        num_threads = config["general_anvio_threads"],
        COGs_dir = config["anvio_COGs_dir"],
        KOs_dir = config["anvio_KOs_dir"]
    log:
        config["logs_dir"] + "make_and_annotate_contigs_db-{genome}.log"
    shell:
        """
        anvi-gen-contigs-database -f {input.fasta} -o {output.contigs_db} -n {params.name} --external-gene-calls {input.gene_calls} -T {params.num_threads} > {log} 2>&1
        anvi-import-functions -c {output.contigs_db} -i {input.gene_functions} > {log} 2>&1

        anvi-run-hmms -T {params.num_threads} -I Bacteria_71 -c {output.contigs_db} > {log} 2>&1

        anvi-scan-trnas -T {params.num_threads} -c {output.contigs_db} > {log} 2>&1

        anvi-run-ncbi-cogs -c {output.contigs_db} --cog-data-dir {params.COGs_dir} -T {params.num_threads} --sensitive > {log} 2>&1

        anvi-run-kegg-kofams -c {output.contigs_db} --kegg-data-dir {params.KOs_dir} -T {params.num_threads} > {log} 2>&1
        """


rule anvi_profile:
    input:
        contigs_db = os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db"),
        bam = os.path.join(config["bam_files_dir"] + "{genome}.bam")
    output:
        profile_db = directory(os.path.join(config["profile_dbs_dir"], "{genome}-profile"))
    params:
        name = "{genome}",
        num_threads = config["general_anvio_threads"]
    log:
        config["logs_dir"] + "anvi_profile-{genome}.log"
    shell:
        """
        anvi-profile -c {input.contigs_db} -i {input.bam} -o {output.profile_db} -S {params.name} \
                     --cluster-contigs --min-contig-length 1000 -T {params.num_threads} > {log} 2>&1
        """


rule gen_external_genomes:
    input:
        contigs_dbs_trigger = expand(os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db"), genome = genome_ID_list),
        profile_dbs_trigger = expand(os.path.join(config["profile_dbs_dir"], "{genome}-profile"), genome = genome_ID_list)
    output:
        "external-genomes.tsv"
    params:
        genomes_file = config["genome_IDs_file"]
    shell:
        """
        printf "name\tcontigs_db_path\n" > {output}

        sed 's/^/contigs-dbs\//' {params.genomes_file} | sed 's/$/-contigs.db/'> paths.tmp

        paste {params.genomes_file} paths.tmp >> {output}

        rm paths.tmp 
        """


rule gen_genomes_storage:
    input:
        "external-genomes.tsv"
    output:
        "our-isolates-GENOMES.db"
    params:
        gene_caller = "NCBI_PGAP"
    log:
        config["logs_dir"] + "gen_genomes_storage.log"
    shell:
        """
        anvi-gen-genomes-storage -e {input} -o {output} --gene-caller {params.gene_caller}
        """


rule anvi_pan:
    input:
        genomes_db = "our-isolates-GENOMES.db"
    output:
        directory(config["anvio_pan_output_dir"])
    params:
        mcl_inflation = config["mcl_inflation"],
        num_threads = config["anvio_pan_threads"],
        name_with_singletons = config["anvio_pan_name"] + "_mcl_" + str(config["mcl_inflation"]),
        name_without_singletons = config["anvio_pan_name"] + "_mcl_" + str(config["mcl_inflation"]) + "_min_2"
    log:
        config["logs_dir"] + "anvi_pan_with_singletons.log"
    shell:
        """
        anvi-pan-genome -g {input.genomes_db} --mcl-inflation {params.mcl_inflation} --min-occurrence 1 \
                        -n {params.name_with_singletons} -o {output} -T {params.num_threads} --sensitive > {log} 2>&1

        printf "\n\n    Now running without singletons:\n\n"

        anvi-pan-genome -g {input.genomes_db} --mcl-inflation {params.mcl_inflation} --min-occurrence 2 \
                        -n {params.name_without_singletons} -o {output} -T {params.num_threads} --sensitive > {log} 2>&1
        """


rule clean_all:
    shell:
        """
        rm -rf {dirs_to_create} "external-genomes.tsv"
        """


