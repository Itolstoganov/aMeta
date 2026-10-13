checkpoint Create_Sample_TaxID_Directories:
    """Create taxid directory

    For a sample, create taxid for each entry in krakenuniq output
    taxID.species. Downstream rules use the taxid directories as
    input, but it is not known beforehand which these are; they are
    determined by the finds in krakenuniq.

    """
    input:
        species="results/KRAKENUNIQ/{sample}/taxID.species",
    output:
        done="results/AUTHENTICATION/{profiler}/{sample}/.extract_taxids_done",
    log:
        "logs/CREATE_SAMPLE_TAXID_DIRECTORIES/{profiler}_{sample}.log",
    threads: 1
    params:
        dir=lambda wildcards: f"results/AUTHENTICATION/{wildcards.profiler}/{wildcards.sample}",
    shell:
        "mkdir -p {params.dir}; "
        "while read taxid; do mkdir -p {params.dir}/$taxid; touch {params.dir}/$taxid/.done; done<{input.species};"
        "touch {output.done}"


rule aggregate:
    """aggregate rule: generate all sample/taxid combinations to
    generate targets.

    The reference id depends on the extraction (Malt_Extract or
    NgsLCA_Extract) having been run, which therefore needs to be triggered
    before the aggregate step. However, beginning from an empty folder
    structure, the extraction output does not exist which means taxid cannot
    be converted to refid. For this reason, the extraction rules are
    checkpoints, which upon completion trigger reevaluation of the DAG.

    """
    input:
        aggregate_extract,
        aggregate_PMD,
        aggregate_plots,
        aggregate_post,
        aggregate_scores,
    output:
        "results/AUTHENTICATION/{profiler}/.{sample}_done",
    log:
        "logs/AGGREGATE/{profiler}_{sample}.log",
    threads: 1
    shell:
        "touch {output}; "


rule Make_Node_List:
    """Generate a list of species names for a taxonomic identifier"""
    input:
        dirdone="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/.done",
    output:
        node_list="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/node_list.txt",
    params:
        tax_db=config["krakenuniq_db"],
    log:
        "logs/MAKE_NODE_LIST/{profiler}_{sample}_{taxid}.log",
    threads: 1
    shell:
        "awk -v var={wildcards.taxid} '{{ if($1==var) print $0 }}' {params.tax_db}/taxDB | cut -f3 > {output.node_list}"


checkpoint Malt_Extract:
    """Convert rma6 output to misc usable formats.

    Downstream rules requires MaltExtract having been run.
    Therefore this rule is a checkpoint that will trigger reevaluation
    of downstream rules. The aggregation is performed by
    aggregate_extract.

    """
    wildcard_constraints:
        profiler="malt",
    input:
        rma6="results/MALT/{sample}.trimmed.rma6",
        node_list="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/node_list.txt",
        ncbi_db_tre=os.path.join(config["ncbi_db"], "ncbi.tre"),
        ncbi_db_map=os.path.join(config["ncbi_db"], "ncbi.map"),
    output:
        extract_log="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/Extract_output/log.txt",
        nodeentries="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/Extract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
    params:
        ncbi_db=lambda wildcards, input: os.path.dirname(input.ncbi_db_tre),
        extract=format_extract_output_directory,
    threads: 4
    log:
        "logs/MALT_EXTRACT/{profiler}_{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    message:
        "Malt_Extract: RUNNING MALT EXTRACT FOR SAMPLE {input.rma6}"
    shell:
        "time MaltExtract -Xmx32G -i {input.rma6} -f def_anc -o {params.extract} -r {params.ncbi_db} --reads --destackingOff --downSampOff --dupRemOff --threads {threads} --matches --minPI 85.0 --maxReadLength 0 --minComp 0.0 --meganSummary -t {input.node_list} -v 2> {log}"


checkpoint NgsLCA_Extract:
    """MaltExtract replacement for aligner + ngsLCA profilers.

    Writes the Extract_output/ tables that score.R and authentic.R read, from
    the aligner's name-sorted BAM and the ngsLCA assignment
    (scripts/ngslca_extract.py). 
    """
    wildcard_constraints:
        profiler=f"({'|'.join(SUPPORTED_ALIGNERS)})",
    input:
        bam="results/ALIGNMENT/{profiler}/{sample}.trimmed.namesorted.bam",
        lca="results/NGSLCA/{profiler}/{sample}.lca",
        node_list="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/node_list.txt",
        ref_fasta="results/REFERENCE_DB/library.project.fna",
        ref_fai="results/REFERENCE_DB/library.project.fna.fai",
        seqid2taxid="results/REFERENCE_DB/seqid2taxid.project.map",
    output:
        extract_log="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/Extract_output/log.txt",
        nodeentries="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/Extract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
    params:
        extract=format_extract_output_directory,
        exe=WORKFLOW_DIR / "scripts/ngslca_extract.py",
        rma6_basename="{sample}.trimmed.rma6",
    threads: 1
    log:
        "logs/NGSLCA_EXTRACT/{profiler}_{sample}_{taxid}.log",
    conda:
        "../envs/ngslca.yaml"
    message:
        "NgsLCA_Extract: BUILDING AUTHENTICATION TABLES FOR SAMPLE {input.bam} TAXID {wildcards.taxid}"
    shell:
        "python {params.exe} --bam {input.bam} --lca {input.lca} "
        "--taxid {wildcards.taxid} --node-list {input.node_list} "
        "--ref-fasta {input.ref_fasta} --seqid2taxid {input.seqid2taxid} "
        "--out-dir {params.extract} --rma6-basename {params.rma6_basename} 2> {log}"


rule Post_Processing:
    input:
        node_list="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/node_list.txt",
    output:
        analysis="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/Extract_output/analysis.RData",
    threads: 4
    params:
        extract=format_extract_output_directory,
    log:
        "logs/POST_PROCESSING/{profiler}_{sample}_{taxid}.log",
    conda:
        "../envs/malt.yaml"
    message:
        "Post_Processing: POSTPROCESSING SAMPLES"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "postprocessing.AMPS.r -m def_anc -r {params.extract} -t {threads} -n {input.node_list} || echo 'postprocessing failed for {wildcards.sample}_{wildcards.taxid}' > {output.analysis}  2> {log}"


rule Samtools_Faidx:
    output:
        fai="{prefix}.{fasta}.fai",
    input:
        fna=ancient("{prefix}.{fasta}"),
    wildcard_constraints:
        fasta="(fna|fasta|fa)",
    message:
        "Samtools_Faidx: INDEXING MALT FASTA DATABASE FOR BREADTH_OF_COVERAGE SEQUENCE RETRIEVAL"
    log:
        "{prefix}.{fasta}.fai_Samtools_Faidx.log",
    threads: 1
    conda:
        "../envs/samtools.yaml"
    envmodules:
        *config["envmodules"]["samtools"],
    shell:
        "samtools faidx {input.fna} 2> {log}"


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


rule Breadth_Of_Coverage:
    # Reference sequences come from the Reference_Subset output, not the full nt FASTA
    input:
        sam=auth_alignment_sam,
        ref_fasta="results/REFERENCE_DB/library.project.fna",
        ref_fasta_fai="results/REFERENCE_DB/library.project.fna.fai",
        nodeentries="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/Extract_output/default/readDist/{sample}.trimmed.rma6_additionalNodeEntries.txt",
    output:
        name_list="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/name_list.txt",
        sorted_bam="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/sorted.bam",
        breadth_of_coverage="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/breadth_of_coverage",
        sam=temporary("results/AUTHENTICATION/{profiler}/{sample}/{taxid}/{taxid}.sam"),
    params:
        ref_id=get_ref_id,
    message:
        "Breadth_Of_Coverage: COMPUTING BREADTH OF COVERAGE, EXTRACTING REFERENCE SEQUENCE FOR VISUALIZING ALIGNMENTS WITH IGV"
    log:
        "logs/BREADTH_OF_COVERAGE/{profiler}_{sample}_{taxid}.log",
    threads: 1
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "echo {params.ref_id} > {output.name_list}; "
        # Keep only records whose RNAME is ref_id.
        # RNAMEs are compared up to the first "|" separator to accomodate the MALT alignment format.
        "zcat -f {input.sam} | awk -F'\\t' -v ref=\"{params.ref_id}\" '"
        "/^@SQ/ {{ split($2, sn, \"|\"); if (sn[1] == \"SN:\" ref) print; next }} "
        "/^@/ {{ print; next }} "
        "{{ split($3, rn, \"|\"); if (rn[1] == ref) print }}"
        "' | uniq > {output.sam}; "
        "samtools view -bS {output.sam} > results/AUTHENTICATION/{wildcards.profiler}/{wildcards.sample}/{wildcards.taxid}/{params.ref_id}.bam; "
        "samtools sort results/AUTHENTICATION/{wildcards.profiler}/{wildcards.sample}/{wildcards.taxid}/{params.ref_id}.bam > {output.sorted_bam}; "
        # Aligner SAMs have no MD tags, which pmdtools needs.
        "if [ \"{wildcards.profiler}\" != \"malt\" ]; then "
        "D=$(dirname {output.sorted_bam}); "
        "samtools faidx {input.ref_fasta} {params.ref_id} > $D/{params.ref_id}.calmd.fasta; "
        "samtools faidx $D/{params.ref_id}.calmd.fasta; "
        "samtools calmd -b {output.sorted_bam} $D/{params.ref_id}.calmd.fasta > $D/sorted.md.bam 2>> {log}; "
        "mv $D/sorted.md.bam {output.sorted_bam}; "
        "fi; "
        "samtools index {output.sorted_bam}; "
        # -g 0x100 also counts secondary alignments.
        "samtools depth -a -g 0x100 {output.sorted_bam} > {output.breadth_of_coverage}; "
        "grep -w -f {output.name_list} {input.ref_fasta_fai} | awk '{{printf(\"%s:1-%s\\n\", $1, $2)}}' > {output.name_list}.regions; "
        "samtools faidx {input.ref_fasta} -r {output.name_list}.regions -o results/AUTHENTICATION/{wildcards.profiler}/{wildcards.sample}/{wildcards.taxid}/{params.ref_id}.fasta"


rule Read_Length_Distribution:
    input:
        bam="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/sorted.bam",
    output:
        distribution="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/read_length.txt",
    message:
        "Read_Length_Distribution: COMPUTING READ LENGTH DISTRIBUTION"
    log:
        "logs/READ_LENGTH_DISTRIBUTION/{profiler}_{sample}_{taxid}.log",
    threads: 1
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "samtools view {input.bam} | awk '{{ print length($10) }}' > {output.distribution}"


rule PMD_scores:
    input:
        bam="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/sorted.bam",
    output:
        scores="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/PMDscores.txt",
    message:
        "PMD_scores: COMPUTING PMD SCORES"
    log:
        "logs/PMD_SCORES/{profiler}_{sample}_{taxid}.log",
    threads: 1
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "(samtools view -h {input.bam} || true) | pmdtools --printDS > {output.scores}"


rule Authentication_Plots:
    input:
        dir="results/AUTHENTICATION/{profiler}/{sample}/{taxid}",
        node_list="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/node_list.txt",
        distribution="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/read_length.txt",
        scores="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/PMDscores.txt",
        breadth_of_coverage="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/breadth_of_coverage",
    output:
        pdf="results/AUTHENTICATION_PLOTS_PDF/{profiler}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf",
        png="results/AUTHENTICATION_PLOTS_PNG/{profiler}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.png",
        pdf_plot="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf",
        png_plot="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.png",
    params:
        exe=WORKFLOW_DIR / "scripts/authentic.R",
    message:
        "Authentication_Plots: MAKING AUTHENTICATION AND VALIDATION PLOTS"
    log:
        "logs/AUTHENTICATION_PLOTS/{profiler}_{sample}_{taxid}.log",
    threads: 1
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "Rscript {params.exe} {wildcards.taxid} {wildcards.sample}.trimmed.rma6 {input.dir}/;"
        "cp {output.pdf_plot} {output.pdf};"
        "cp {output.png_plot} {output.png};"


rule Deamination:
    input:
        bam="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/sorted.bam",
    output:
        tmp="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/PMD_temp.txt",
        pmd="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/PMD_plot.frag.pdf",
    message:
        "Deamination: INFERRING DEAMINATION PATTERN WITH PMDTOOLS"
    log:
        "logs/DEAMINATION/{profiler}_{sample}_{taxid}.log",
    threads: 1
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "(samtools view {input.bam} || true) | pmdtools --platypus --number 2000000 > {output.tmp}; "
        "cd results/AUTHENTICATION/{wildcards.profiler}/{wildcards.sample}/{wildcards.taxid}; "
        "R CMD BATCH $(which plotPMD); "


rule Authentication_Score:
    input:
        extract_log="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/Extract_output/log.txt",
        name_list="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/name_list.txt",
        scores="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/PMDscores.txt",
        breadth_of_coverage="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/breadth_of_coverage",
        read_length="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/read_length.txt",
    output:
        scores="results/AUTHENTICATION/{profiler}/{sample}/{taxid}/authentication_scores.txt",
    message:
        "Authentication_Score: COMPUTING AUTHENTICATION SCORES"
    params:
        exe=WORKFLOW_DIR / "scripts/score.R",
        # score.R only uses the rma6 file name, not the contents.
        rma6_basename="{sample}.trimmed.rma6",
    log:
        "logs/AUTHENTICATION_SCORE/{profiler}_{sample}_{taxid}.log",
    threads: 1
    conda:
        "../envs/malt.yaml"
    envmodules:
        *config["envmodules"]["malt"],
    shell:
        "Rscript {params.exe} {params.rma6_basename} $(dirname {input.extract_log}) {input.name_list} $(dirname {input.name_list}) {input.scores} &> {log};"
