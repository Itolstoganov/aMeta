BOWTIE2_INDEX_EXT = [".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"]


def bowtie2_index(fasta):
    """bowtie2-build output for `fasta`, written next to it with the FASTA path as prefix."""
    return [fasta + ext for ext in BOWTIE2_INDEX_EXT]


rule Bowtie2_Index:
    """Index a FASTA for bowtie2: bowtie2_db for the pathogen screen, and the
    project reference for the bowtie2 + ngsLCA profiler."""
    output:
        bowtie2_index("{fasta}"),
    input:
        ref=lambda wildcards: (
            ancient(wildcards.fasta)
            if wildcards.fasta == config["bowtie2_db"]
            else wildcards.fasta
        ),
    wildcard_constraints:
        fasta=r".+(?<!\.rev)",
    conda:
        "../envs/bowtie2.yaml"
    envmodules:
        *config["envmodules"]["bowtie2"],
    # single-threaded for bowtie2_db to prioritize memory over runtime
    threads: lambda wildcards: 1 if wildcards.fasta == config["bowtie2_db"] else 20
    log:
        "{fasta}_BOWTIE2_BUILD.log",
    benchmark:
        "{fasta}_BOWTIE2_BUILD.benchmark.txt"
    message:
        "Bowtie2_Index: INDEXING {input.ref} WITH BOWTIE2"
    shell:
        "bowtie2-build --large-index --threads {threads} {input.ref} {input.ref} > {log} 2>&1"


rule Bowtie2_Alignment:
    output:
        bam="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam",
        bai="results/BOWTIE2/{sample}/AlignedToBowtie2DB.bam.bai",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        db=bowtie2_index(config["bowtie2_db"]),
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


def strobealign_binary(wildcards):
    return "resources/bin/strobealign/" + config["aligners"]["strobealign"]["ref"]


rule Strobealign_Index:
    """Write the strobemer index of the project reference once, for all samples.

    The index depends on args and on the read-length profile, so both rules
    pass the same args and read_length; --use-index refuses a mismatched index.
    """
    output:
        "results/REFERENCE_DB/library.project.fna.sti",
    input:
        ref="results/REFERENCE_DB/library.project.fna",
        binary=strobealign_binary,
    params:
        args=strobealign_setting("args"),
        read_length=strobealign_setting("read_length"),
    threads: 20
    log:
        "logs/STROBEALIGN_INDEX/STROBEALIGN_INDEX.log",
    benchmark:
        "benchmarks/STROBEALIGN_INDEX/STROBEALIGN_INDEX.benchmark.txt"
    message:
        "Strobealign_Index: INDEXING THE PROJECT REFERENCE WITH STROBEALIGN"
    shell:
        "{input.binary} {params.args} -r {params.read_length} -t {threads} "
        "--create-index {input.ref} > {log} 2>&1"


rule Strobealign_Alignment:
    """Align trimmed reads with strobealign, dropping short local alignments
    (scripts/filter_alignments.py)."""
    output:
        sam="results/ALIGNMENT/strobealign/{sample}.trimmed.sam.gz",
        namesorted_bam="results/ALIGNMENT/strobealign/{sample}.trimmed.namesorted.bam",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        ref="results/REFERENCE_DB/library.project.fna",
        index="results/REFERENCE_DB/library.project.fna.sti",
        binary=strobealign_binary,
    params:
        sam="results/ALIGNMENT/strobealign/{sample}.trimmed.sam",
        args=strobealign_setting("args"),
        read_length=strobealign_setting("read_length"),
        min_aligned_bp=config["alignment_filter"]["min_aligned_bp"],
        min_aligned_frac=config["alignment_filter"]["min_aligned_frac"],
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
        "{input.binary} {params.args} -r {params.read_length} -t {threads} "
        "--use-index {input.ref} {input.fastq} 2> {log} "
        "| python {params.filter_script} {params.min_aligned_bp} {params.min_aligned_frac} "
        "2>> {log} > {params.sam}; "
        "pigz -c {params.sam} > {output.sam}; "
        "samtools sort -n -@ {threads} -o {output.namesorted_bam} {params.sam}; "
        "rm {params.sam}"


def bowtie2_setting(key, default=None):
    return lambda wildcards: config["aligners"]["bowtie2"].get(key, default)


rule Bowtie2_Profiler_Alignment:
    """Align trimmed reads with bowtie2 end-to-end, reporting up to -k hits per read.

    End-to-end mode never soft-clips, so terminal damage stays in the alignment.
    Hits more than max_score_drop below the read's best AS are dropped, as
    strobealign's --st does for its secondary alignments.
    """
    output:
        sam="results/ALIGNMENT/bowtie2/{sample}.trimmed.sam.gz",
        namesorted_bam="results/ALIGNMENT/bowtie2/{sample}.trimmed.namesorted.bam",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        ref="results/REFERENCE_DB/library.project.fna",
        index=bowtie2_index("results/REFERENCE_DB/library.project.fna"),
    params:
        sam="results/ALIGNMENT/bowtie2/{sample}.trimmed.sam",
        args=bowtie2_setting("args"),
        min_aligned_bp=config["alignment_filter"]["min_aligned_bp"],
        min_aligned_frac=config["alignment_filter"]["min_aligned_frac"],
        max_score_drop=bowtie2_setting("max_score_drop", ""),
        filter_script=WORKFLOW_DIR / "scripts/filter_alignments.py",
    threads: 20
    log:
        "logs/ALIGNMENT/bowtie2/{sample}.log",
    conda:
        "../envs/aligner.yaml"
    envmodules:
        *config["envmodules"]["bowtie2"],
    benchmark:
        "benchmarks/ALIGNMENT/bowtie2/{sample}.benchmark.txt"
    message:
        "Bowtie2_Profiler_Alignment: ALIGNING SAMPLE {input.fastq} WITH BOWTIE2"
    shell:
        "bowtie2 {params.args} --no-unal -p {threads} -x {input.ref} "
        "-U {input.fastq} 2> {log} "
        "| python {params.filter_script} {params.min_aligned_bp} {params.min_aligned_frac} "
        "{params.max_score_drop} 2>> {log} > {params.sam}; "
        "pigz -c {params.sam} > {output.sam}; "
        "samtools sort -n -@ {threads} -o {output.namesorted_bam} {params.sam}; "
        "rm {params.sam}"


def bwa_aln_setting(key, default=None):
    return lambda wildcards: config["aligners"]["bwa_aln"].get(key, default)


BWA_INDEX_EXT = [".amb", ".ann", ".bwt", ".pac", ".sa"]


rule Bwa_Index:
    """Index the project reference for the bwa aln + ngsLCA profiler."""
    output:
        ["results/REFERENCE_DB/library.project.fna" + ext for ext in BWA_INDEX_EXT],
    input:
        ref="results/REFERENCE_DB/library.project.fna",
    threads: 1
    log:
        "logs/BWA_INDEX/BWA_INDEX.log",
    conda:
        "../envs/aligner.yaml"
    benchmark:
        "benchmarks/BWA_INDEX/BWA_INDEX.benchmark.txt"
    message:
        "Bwa_Index: INDEXING THE PROJECT REFERENCE WITH BWA"
    shell:
        "bwa index {input.ref} > {log} 2>&1"


rule Bwa_Aln_Profiler_Alignment:
    """Align trimmed reads with bwa aln + samse, writing each XA hit as a
    secondary record (scripts/bwa_xa_to_secondary.py).

    bwa aln aligns end to end, so terminal damage stays in the alignment.
    """
    output:
        sam="results/ALIGNMENT/bwa_aln/{sample}.trimmed.sam.gz",
        namesorted_bam="results/ALIGNMENT/bwa_aln/{sample}.trimmed.namesorted.bam",
    input:
        fastq="results/CUTADAPT_ADAPTER_TRIMMING/{sample}.trimmed.fastq.gz",
        ref="results/REFERENCE_DB/library.project.fna",
        index=rules.Bwa_Index.output,
    params:
        sai="results/ALIGNMENT/bwa_aln/{sample}.sai",
        sam="results/ALIGNMENT/bwa_aln/{sample}.trimmed.sam",
        args=bwa_aln_setting("args"),
        max_hits=bwa_aln_setting("max_hits"),
        min_aligned_bp=config["alignment_filter"]["min_aligned_bp"],
        min_aligned_frac=config["alignment_filter"]["min_aligned_frac"],
        xa_script=WORKFLOW_DIR / "scripts/bwa_xa_to_secondary.py",
        filter_script=WORKFLOW_DIR / "scripts/filter_alignments.py",
    threads: 20
    log:
        "logs/ALIGNMENT/bwa_aln/{sample}.log",
    conda:
        "../envs/aligner.yaml"
    benchmark:
        "benchmarks/ALIGNMENT/bwa_aln/{sample}.benchmark.txt"
    message:
        "Bwa_Aln_Profiler_Alignment: ALIGNING SAMPLE {input.fastq} WITH BWA ALN"
    shell:
        "bwa aln {params.args} -t {threads} {input.ref} {input.fastq} > {params.sai} 2> {log}; "
        "bwa samse -n {params.max_hits} {input.ref} {params.sai} {input.fastq} 2>> {log} "
        "| python {params.xa_script} {params.max_hits} 2>> {log} "
        "| python {params.filter_script} {params.min_aligned_bp} {params.min_aligned_frac} "
        "2>> {log} > {params.sam}; "
        "pigz -c {params.sam} > {output.sam}; "
        "samtools sort -n -@ {threads} -o {output.namesorted_bam} {params.sam}; "
        "rm {params.sai} {params.sam}"
