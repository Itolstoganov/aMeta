rule Malt_Acc2Taxa:
    """Convert NCBI's accession2taxid into the two-column form malt-build reads.

    malt-build's -a2t parser takes field 1 of each line as the accession and
    field 2 as the taxon id. NCBI's four-column nucl_gb.accession2taxid is not parsed correctly
    """
    output:
        acc2taxa="results/MALT_DB/acc2taxa.tsv",
    input:
        accession2taxid=config["malt_accession2taxid"],
    log:
        "logs/BUILD_MALT_DB/MALT_ACC2TAXA.log",
    threads: 1
    message:
        "Malt_Acc2Taxa: CONVERTING {input.accession2taxid} TO malt-build -a2t FORMAT"
    shell:
        # -f passes an uncompressed file straight through, so this takes the
        # NCBI file either as shipped (.gz) or already decompressed. The field 3
        # test drops the header line and any malformed row.
        "gzip -cdf {input.accession2taxid} | "
        "awk -F'\t' '$3 ~ /^[0-9]+$/ {{ print $1 \"\t\" $3 }}' > {output.acc2taxa} 2> {log}"


rule Build_Malt_DB:
    output:
        db=temp(directory("results/MALT_DB/maltDB.dat")),
    input:
        # The FASTA subset is prepared by the shared Reference_Subset rule
        # (workflow/rules/align.smk); this rule only runs malt-build on it.
        project_fasta="results/REFERENCE_DB/library.project.fna",
        # two-column form of config["malt_accession2taxid"]; see Malt_Acc2Taxa
        acc2taxa="results/MALT_DB/acc2taxa.tsv",
    params:
        max_heap=get_malt_heap_option(),
        # malt-build labels a reference by parsing organism names out of its
        # FASTA description as well as looking the accession up in -a2t, and
        # takes the LCA of every name it matches -- so a description naming a
        # second organism (a homolog it is "similar to", a cloning host, a lab
        # in a clone name) drags the label to a coarser or unrelated node. Off,
        # labels come from the accession alone, which is what ngsLCA's -acc2tax
        # does, making the two backends' taxon assignments comparable.
        parse_taxon_names=config.get("malt_parse_taxon_names", True),
    threads: 20
    log:
        "logs/BUILD_MALT_DB/BUILD_MALT_DB.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    benchmark:
        "benchmarks/BUILD_MALT_DB/BUILD_MALT_DB.benchmark.txt"
    message:
        "Build_Malt_DB: BUILDING MALT DATABASE USING SPECIES DETECTED BY KRAKENUNIQ"
    script:
        "../scripts/malt-build.py"


rule Malt:
    output:
        rma6="results/MALT/{sample}.trimmed.rma6",
        sam="results/MALT/{sample}.trimmed.sam.gz",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        db="results/MALT_DB/maltDB.dat",
    params:
        gunzipped_sam="results/MALT/{sample}.trimmed.sam",
        max_heap=get_malt_heap_option(),
    threads: 20
    log:
        "logs/MALT/{sample}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    benchmark:
        "benchmarks/MALT/{sample}.benchmark.txt"
    message:
        "Malt: RUNNING MALT ALIGNMENTS FOR SAMPLE {input.fastq}"
    shell:
        "unset DISPLAY; malt-run {params.max_heap} -at SemiGlobal -m BlastN -i {input.fastq} -o {output.rma6} -a {params.gunzipped_sam} -t {threads} -d {input.db} -sup 1 -mq 100 -top 1 -mpi 85.0 -id 85.0 -v &> {log}"


rule Malt_QuantifyAbundance:
    output:
        out_dir=directory("results/MALT_QUANTIFY_ABUNDANCE/{sample}"),
        counts="results/MALT_QUANTIFY_ABUNDANCE/{sample}/sam_counts.txt",
    input:
        sam="results/MALT/{sample}.trimmed.sam.gz",
    params:
        unique_taxids="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt",
        exe=WORKFLOW_DIR / "scripts/malt_quantify_abundance.py",
    log:
        "logs/MALT_QUANTIFY_ABUNDANCE/{sample}.log",
    threads: 1
    benchmark:
        "benchmarks/MALT_QUANTIFY_ABUNDANCE/{sample}.benchmark.txt"
    message:
        "Malt_QuantifyAbundance: QUANTIFYING MICROBIAL ABUNDANCE USING MALT SAM-ALIGNMENTS FOR SAMPLE {input.sam}"
    shell:
        "mkdir -p {output.out_dir}; "
        "{params.exe} {input.sam} {params.unique_taxids} > {output.counts} 2> {log}"


rule Malt_AbundanceMatrix_Sam:
    output:
        out_dir=directory("results/MALT_ABUNDANCE_MATRIX_SAM"),
        abundance_matrix="results/MALT_ABUNDANCE_MATRIX_SAM/malt_abundance_matrix_sam.txt",
    input:
        sam_counts=expand(
            "results/MALT_QUANTIFY_ABUNDANCE/{sample}/sam_counts.txt", sample=SAMPLES
        ),
    log:
        "logs/MALT_ABUNDANCE_MATRIX_SAM/MALT_ABUNDANCE_MATRIX_SAM.log",
    threads: 1
    params:
        exe=WORKFLOW_DIR / "scripts/malt_abundance_matrix.R",
    conda:
        "../envs/r.yaml"
    envmodules:
        *config["envmodules"]["r"],
    benchmark:
        "benchmarks/MALT_ABUNDANCE_MATRIX_SAM/MALT_ABUNDANCE_MATRIX_SAM.benchmark.txt"
    message:
        "Malt_AbundanceMatrix_Sam: COMPUTING MALT MICROBIAL ABUNDANCE MATRIX FROM SAM-FILES"
    shell:
        "Rscript {params.exe} results/MALT_QUANTIFY_ABUNDANCE {output.out_dir} &> {log}"


rule Malt_AbundanceMatrix_Rma6:
    output:
        out_dir=directory("results/MALT_ABUNDANCE_MATRIX_RMA6"),
        abundance_matrix="results/MALT_ABUNDANCE_MATRIX_RMA6/malt_abundance_matrix_rma6.txt",
    input:
        rma6=expand("results/MALT/{sample}.trimmed.rma6", sample=SAMPLES),
    params:
        exe=WORKFLOW_DIR / "scripts/rma-tabuliser",
    log:
        "logs/MALT_ABUNDANCE_MATRIX_RMA6/MALT_ABUNDANCE_MATRIX_RMA6.log",
    threads: 1
    envmodules:
        *config["envmodules"]["malt"],
    conda:
        "../envs/malt.yaml"
    benchmark:
        "benchmarks/MALT_ABUNDANCE_MATRIX_RMA6/MALT_ABUNDANCE_MATRIX_RMA6.benchmark.txt"
    message:
        "Malt_AbundanceMatrix_Rma6: COMPUTING MALT MICROBIAL ABUNDANCE MATRIX FROM RMA6-FILES"
    shell:
        "{params.exe} -d $(dirname {input.rma6}) -r 'S' &> {log}; "
        "mv results/MALT/count_table.tsv {output.out_dir}; "
        "mv {output.out_dir}/count_table.tsv {output.abundance_matrix}"


if pv.parse(snakemake_version) >= pv.parse("8.0.0"):
    rule NCBIMapTre:
        """Download ncbi.map and ncbi.tre from https://github.com/husonlab/megan-ce/tree/master/src/megan/resources/files"""
        output:
            tre=os.path.join(config["ncbi_db"], "ncbi.tre"),
            map=os.path.join(config["ncbi_db"], "ncbi.map"),
        input:
            tre=storage(
                "https://github.com/husonlab/megan-ce/raw/master/src/megan/resources/files/ncbi.tre",
            ),
            map=storage(
                "https://github.com/husonlab/megan-ce/raw/master/src/megan/resources/files/ncbi.map",
            ),
        log:
            "logs/NCBI/ncbi.log",
        threads: 1
        shell:
            "mv {input.tre} {output.tre};"
            "mv {input.map} {output.map};"
else:
    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule NCBIMapTre:
        """Download ncbi.map and ncbi.tre from https://github.com/husonlab/megan-ce/tree/master/src/megan/resources/files"""
        output:
            tre=os.path.join(config["ncbi_db"], "ncbi.tre"),
            map=os.path.join(config["ncbi_db"], "ncbi.map"),
        input:
            tre=HTTP.remote(
                "github.com/husonlab/megan-ce/raw/master/src/megan/resources/files/ncbi.tre",
                keep_local=True,
            ),
            map=HTTP.remote(
                "github.com/husonlab/megan-ce/raw/master/src/megan/resources/files/ncbi.map",
                keep_local=True,
            ),
        log:
            "logs/NCBI/ncbi.log",
        threads: 1
        shell:
            "mv {input.tre} {output.tre};"
            "mv {input.map} {output.map};"
