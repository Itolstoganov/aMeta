## ngsLCA taxonomy (aligner + ngsLCA profilers)
##
## ngsLCA assigns each read a lowest-common-ancestor taxid from the aligner's
## name-sorted BAM. Outputs are per profiler, under results/NGSLCA/<profiler>/.
## names.dmp/nodes.dmp are derived from the KrakenUniq taxDB; -acc2tax is
## config["malt_accession2taxid"] (NCBI accession2taxid).


rule KrakenUniq_TaxDB_To_Dmp:
    """Convert the KrakenUniq taxDB into NCBI names.dmp/nodes.dmp for ngsLCA."""
    output:
        names="results/NGSLCA_DB/names.dmp",
        nodes="results/NGSLCA_DB/nodes.dmp",
    input:
        taxdb=os.path.join(config["krakenuniq_db"], "taxDB"),
    params:
        exe=WORKFLOW_DIR / "scripts/taxdb_to_dmp.py",
    threads: 1
    log:
        "logs/NGSLCA_DB/taxdb_to_dmp.log",
    conda:
        "../envs/ngslca.yaml"
    message:
        "KrakenUniq_TaxDB_To_Dmp: CONVERTING KRAKENUNIQ taxDB TO names.dmp/nodes.dmp"
    shell:
        "python {params.exe} {input.taxdb} {output.names} {output.nodes} 2> {log}"


rule NgsLCA:
    """Assign a lowest common ancestor taxid to every read with ngsLCA."""
    output:
        lca="results/NGSLCA/{profiler}/{sample}.lca",
    input:
        bam="results/ALIGNMENT/{profiler}/{sample}.trimmed.namesorted.bam",
        names="results/NGSLCA_DB/names.dmp",
        nodes="results/NGSLCA_DB/nodes.dmp",
        acc2tax=config["malt_accession2taxid"],
    params:
        prefix="results/NGSLCA/{profiler}/{sample}",
        args=lambda wildcards: config["ngslca_args"],
    threads: 4
    log:
        "logs/NGSLCA/{profiler}/{sample}.log",
    conda:
        "../envs/ngslca.yaml"
    benchmark:
        "benchmarks/NGSLCA/{profiler}/{sample}.benchmark.txt"
    message:
        "NgsLCA: ASSIGNING TAXONOMY WITH ngsLCA FOR SAMPLE {input.bam}"
    shell:
        "ngsLCA {params.args} -names {input.names} -nodes {input.nodes} "
        "-acc2tax {input.acc2tax} -bam {input.bam} -outnames {params.prefix} 2> {log}"


rule NgsLCA_AbundanceMatrix:
    """Count reads per LCA taxid across samples into an abundance matrix."""
    output:
        abundance_matrix="results/NGSLCA_ABUNDANCE_MATRIX/{profiler}/ngslca_abundance_matrix.txt",
    input:
        lca=expand("results/NGSLCA/{{profiler}}/{sample}.lca", sample=SAMPLES),
    params:
        exe=WORKFLOW_DIR / "scripts/ngslca_abundance_matrix.py",
        pairs=lambda wildcards, input: " ".join(
            f"{s} {f}" for s, f in zip(SAMPLES, input.lca)
        ),
    threads: 1
    log:
        "logs/NGSLCA_ABUNDANCE_MATRIX/{profiler}.log",
    conda:
        "../envs/ngslca.yaml"
    message:
        "NgsLCA_AbundanceMatrix: BUILDING ngsLCA ABUNDANCE MATRIX FOR PROFILER {wildcards.profiler}"
    shell:
        "python {params.exe} {output.abundance_matrix} {params.pairs} 2> {log}"
