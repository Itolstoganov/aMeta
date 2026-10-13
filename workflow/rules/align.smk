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


rule Reference_Subset:
    """Subset the nt FASTA to the KrakenUniq-detected species."""
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



## ---------------------------------------------------------------------------
## Read aligner rules for aligner + ngsLCA profilers.
##
## Settings come from config["aligners"][<aligner>]. To add an aligner, add its
## rules here and its name to SUPPORTED_ALIGNERS in common.smk.
## ---------------------------------------------------------------------------

# Settings are read lazily, so runs without strobealign need no aligners.strobealign.
def strobealign_setting(key):
    return lambda wildcards: config["aligners"]["strobealign"][key]


rule Compile_Strobealign:
    """Build the Rust aDNA strobealign from a GitHub archive of the requested ref."""
    output:
        binary="resources/bin/strobealign/{ref}",
    params:
        repo=strobealign_setting("repo"),
        ref="{ref}",
    threads: 4
    log:
        "logs/COMPILE_STROBEALIGN/{ref}.log",
    conda:
        "../envs/strobealign-build.yaml"
    benchmark:
        "benchmarks/COMPILE_STROBEALIGN/{ref}.benchmark.txt"
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
    """Align trimmed reads with strobealign, dropping short local alignments
    (scripts/filter_alignments.py)."""
    output:
        sam="results/ALIGNMENT/strobealign/{sample}.trimmed.sam.gz",
        namesorted_bam="results/ALIGNMENT/strobealign/{sample}.trimmed.namesorted.bam",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        ref="results/REFERENCE_DB/library.project.fna",
        binary=lambda wildcards: "resources/bin/strobealign/"
        + config["aligners"]["strobealign"]["ref"],
    params:
        sam="results/ALIGNMENT/strobealign/{sample}.trimmed.sam",
        args=strobealign_setting("args"),
        min_aligned_bp=strobealign_setting("min_aligned_bp"),
        min_aligned_frac=strobealign_setting("min_aligned_frac"),
        filter_script=WORKFLOW_DIR / "scripts/filter_alignments.py",
    threads: 20
    log:
        "logs/ALIGNMENT/strobealign/{sample}.log",
    conda:
        "../envs/aligner.yaml"
    envmodules:
        *config["envmodules"]["samtools"],
    benchmark:
        "benchmarks/ALIGNMENT/strobealign/{sample}.benchmark.txt"
    message:
        "Strobealign_Alignment: ALIGNING SAMPLE {input.fastq} WITH STROBEALIGN"
    shell:
        "{input.binary} {params.args} -t {threads} "
        "{input.ref} {input.fastq} 2> {log} "
        "| python {params.filter_script} {params.min_aligned_bp} {params.min_aligned_frac} "
        "2>> {log} > {params.sam}; "
        "pigz -c {params.sam} > {output.sam}; "
        "samtools sort -n -@ {threads} -o {output.namesorted_bam} {params.sam}; "
        "rm {params.sam}"
