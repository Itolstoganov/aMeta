rule Bowtie2_Index:
    output:
        expand(
            f"{config['bowtie2_db']}{{ext}}",
            ext=[
                ".1.bt2l",
                ".2.bt2l",
                ".3.bt2l",
                ".4.bt2l",
                ".rev.1.bt2l",
                ".rev.2.bt2l",
            ],
        ),
    input:
        ref=ancient(config["bowtie2_db"]),
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        *config["envmodules"]["bowtie2"],
    threads: 1
    log:
        f"{config['bowtie2_db']}_BOWTIE2_BUILD.log",
    shell:
        "bowtie2-build --large-index --threads {threads} {input.ref} {input.ref} > {log} 2>&1"


rule Bowtie2_Alignment:
    output:
        bam="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam",
        bai="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam.bai",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        db=rules.Bowtie2_Index.output,
    params:
        BOWTIE2_DB=lambda wildcards, input: config["bowtie2_db"],
    threads: 10
    log:
        "logs/BOWTIE2/{sample}.log",
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        *config["envmodules"]["bowtie2"],
    benchmark:
        "benchmarks/BOWTIE2/{sample}.benchmark.txt"
    message:
        "Bowtie2_Alignment: ALIGNING SAMPLE {input.fastq} WITH BOWTIE2"
    shell:
        """bowtie2 --large-index -x {params.BOWTIE2_DB} --end-to-end --threads {threads} --very-sensitive -U {input.fastq} > $(dirname {output.bam})/AlignedToBowtie2DB.sam 2> {log};"""
        """grep @ $(dirname {output.bam})/AlignedToBowtie2DB.sam | awk '!seen[$2]++' > $(dirname {output.bam})/header_nodups.txt;"""
        """grep -v '^@' $(dirname {output.bam})/AlignedToBowtie2DB.sam > $(dirname {output.bam})/AlignedToBowtie2DB.noheader.sam;"""
        """cat $(dirname {output.bam})/header_nodups.txt $(dirname {output.bam})/AlignedToBowtie2DB.noheader.sam > $(dirname {output.bam})/AlignedToBowtie2DB.nodups.sam;"""
        """samtools view -bS -q 1 -h -@ {threads} $(dirname {output.bam})/AlignedToBowtie2DB.nodups.sam | samtools sort - -@ {threads} -o {output.bam};"""
        """samtools index {output.bam};"""
        """rm $(dirname {output.bam})/header_nodups.txt;"""
        """rm $(dirname {output.bam})/AlignedToBowtie2DB.noheader.sam;"""
        """rm $(dirname {output.bam})/AlignedToBowtie2DB.nodups.sam;"""
        """rm $(dirname {output.bam})/AlignedToBowtie2DB.sam;"""


## ---------------------------------------------------------------------------
## Pluggable aligner framework (alternative to MALT alignment)
##
## When any aligner+ngsLCA profiler is configured (NGSLCA_PROFILERS non-empty), a
## read aligner replaces MALT's alignment step for that profiler. The taxonomy is
## then assigned by ngsLCA (see ngslca.smk).
## Everything downstream of the aligner consumes two canonical, aligner-agnostic
## outputs, so the rest of the pipeline is unaware of which aligner produced them:
##
##   results/ALIGNMENT/{sample}.trimmed.sam.gz          SAM, same accession|tax|taxid
##                                                      RNAME space as the MALT SAM
##   results/ALIGNMENT/{sample}.trimmed.namesorted.bam  name-sorted BAM for ngsLCA
##
## Adding another aligner = add an `elif ALIGNER == "..."` block below whose rule
## produces those two paths; no downstream rule needs to change.
## ---------------------------------------------------------------------------


rule Reference_Subset:
    """Subset the nt FASTA to the KrakenUniq-detected species (aligner-agnostic).

    Shared reference preparation used by BOTH backends: Build_Malt_DB feeds this
    subset to malt-build, and the read aligners index it directly. Defined
    unconditionally so either profiler can depend on it. (This is the FASTA
    subsetting that malt-build.py used to do inline.)
    """
    output:
        seqid2taxid_project="results/REFERENCE_DB/seqid2taxid.project.map",
        seqids_project="results/REFERENCE_DB/seqids.project",
        project_headers="results/REFERENCE_DB/project.headers",
        project_fasta="results/REFERENCE_DB/library.project.fna",
    input:
        unique_taxids="results/KRAKENUNIQ_ABUNDANCE_MATRIX/unique_species_taxid_list.txt",
    params:
        seqid2taxid=config["malt_seqid2taxid_db"],
        nt_fasta=config["malt_nt_fasta"],
    threads: 1
    log:
        "logs/REFERENCE_SUBSET/REFERENCE_SUBSET.log",
    conda:
        "../envs/aligner.yaml"
    envmodules:
        *config["envmodules"]["samtools"],
    benchmark:
        "benchmarks/REFERENCE_SUBSET/REFERENCE_SUBSET.benchmark.txt"
    message:
        "Reference_Subset: SUBSETTING NT FASTA TO KRAKENUNIQ-DETECTED SPECIES"
    shell:
        "grep -wFf {input.unique_taxids} {params.seqid2taxid} > {output.seqid2taxid_project}; "
        "cut -f1 {output.seqid2taxid_project} > {output.seqids_project}; "
        "grep -Ff {output.seqids_project} {params.nt_fasta} | sed 's/>//g' > {output.project_headers}; "
        "seqtk subseq {params.nt_fasta} {output.project_headers} > {output.project_fasta} 2> {log}"


rule Project_Fasta_Faidx:
    """Index the project FASTA; Samtools_Faidx's ancient() input would keep the previous run's stale .fai."""
    output:
        fai="results/REFERENCE_DB/library.project.fna.fai",
    input:
        fna="results/REFERENCE_DB/library.project.fna",
    threads: 1
    log:
        "logs/REFERENCE_SUBSET/PROJECT_FASTA_FAIDX.log",
    conda:
        "../envs/samtools.yaml"
    envmodules:
        *config["envmodules"]["samtools"],
    message:
        "Project_Fasta_Faidx: INDEXING THE PROJECT FASTA SUBSET"
    shell:
        "samtools faidx {input.fna} 2> {log}"


ruleorder: Project_Fasta_Faidx > Samtools_Faidx


if NGSLCA_PROFILERS:

    _aligner = config.get("aligner", {})
    ALIGNER = _aligner.get("name", "strobealign")

    if ALIGNER == "strobealign":

        STROBEALIGN_REPO = _aligner.get("repo", "https://github.com/ksahlin/strobealign")
        STROBEALIGN_REF = _aligner.get("ref", "cc24cbd434f4ea7dd3c7fd344e3aaaf774cf4c13")
        STROBEALIGN_ARGS = _aligner.get("args", "--adna -k 17 -s 13 --ry-len 8 -M 500 --ssw -L 32")
        STROBEALIGN_MAX_SECONDARY = _aligner.get("max_secondary", 200)
        STROBEALIGN_BIN = f"resources/bin/strobealign/{STROBEALIGN_REF}"

        rule Compile_Strobealign:
            """Build the (Rust) aDNA strobealign from source at the configured ref.

            The aDNA strobealign (--ancient-dna) is a Rust implementation not available
            on conda, so it is compiled from a configurable commit/branch, mirroring the
            top-level sim benchmark (workflow/rules/alignment.smk::compile_strobealign):
            download the GitHub archive, `cargo build`, install the debug binary. Cached
            under resources/ and reused across all samples.
            """
            output:
                binary=STROBEALIGN_BIN,
            params:
                repo=STROBEALIGN_REPO,
                ref=STROBEALIGN_REF,
            threads: 4
            log:
                "logs/COMPILE_STROBEALIGN/compile_strobealign.log",
            conda:
                "../envs/strobealign-build.yaml"
            benchmark:
                "benchmarks/COMPILE_STROBEALIGN/compile_strobealign.benchmark.txt"
            message:
                "Compile_Strobealign: BUILDING STROBEALIGN {params.ref} FROM SOURCE (cargo)"
            shell:
                "mkdir -p $(dirname {output.binary}) $(dirname {log}); "
                "OUT=$(readlink -f {output.binary}); "
                "LOG=$(readlink -f {log}); "
                "tmp=$(mktemp -d); "
                "trap 'rm -rf $tmp' EXIT; "
                "wget -q {params.repo}/archive/{params.ref}.zip -O $tmp/strobealign.zip 2> $LOG; "
                "unzip -q $tmp/strobealign.zip -d $tmp; "
                "cd $tmp/strobealign-*; "
                "cargo build -j {threads} >> $LOG 2>&1; "
                "mv target/debug/strobealign $OUT"

        rule Strobealign_Alignment:
            """Align trimmed reads with the compiled aDNA strobealign.

            `-N` retains secondary alignments so ngsLCA can resolve the lowest common
            ancestor across every reference a read maps to. Writes the canonical SAM
            (accession|tax|taxid RNAME space) and a name-sorted BAM for ngsLCA.
            """
            output:
                sam="results/ALIGNMENT/{sample}.trimmed.sam.gz",
                namesorted_bam="results/ALIGNMENT/{sample}.trimmed.namesorted.bam",
            input:
                fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
                ref="results/REFERENCE_DB/library.project.fna",
                binary=STROBEALIGN_BIN,
            params:
                sam="results/ALIGNMENT/{sample}.trimmed.sam",
                args=STROBEALIGN_ARGS,
                max_secondary=STROBEALIGN_MAX_SECONDARY,
            threads: 20
            log:
                "logs/ALIGNMENT/{sample}.log",
            conda:
                "../envs/aligner.yaml"
            envmodules:
                *config["envmodules"]["samtools"],
            benchmark:
                "benchmarks/ALIGNMENT/{sample}.benchmark.txt"
            message:
                "Strobealign_Alignment: ALIGNING SAMPLE {input.fastq} WITH STROBEALIGN"
            shell:
                "{input.binary} {params.args} -N {params.max_secondary} -t {threads} "
                "{input.ref} {input.fastq} 2> {log} > {params.sam}; "
                "pigz -c {params.sam} > {output.sam}; "
                "samtools sort -n -@ {threads} -o {output.namesorted_bam} {params.sam}; "
                "rm {params.sam}"

    else:
        raise WorkflowError(
            f"Unknown aligner '{ALIGNER}'. Supported: strobealign. "
            "Add an `elif ALIGNER == ...` block in workflow/rules/align.smk."
        )
