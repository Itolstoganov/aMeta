## ngsLCA taxonomy rules (used when an aligner+ngsLCA profiler is configured).
##
## ngsLCA (https://github.com/miwipe/ngsLCA) assigns each read a lowest-common-
## ancestor taxid from a name-sorted BAM (produced by the aligner, see align.smk)
## plus an NCBI taxonomy. It replaces MALT's LCA step; the per-read `.lca` output
## feeds both the abundance matrix (here) and the authentication tables (see
## ngslca_extract.py, wired in authentic.smk).
##
## Taxonomy inputs reuse the existing aMeta databases:
##   - names.dmp / nodes.dmp are derived from the KrakenUniq taxDB (which ships no
##     dmp files) by KrakenUniq_TaxDB_To_Dmp.
##   - -acc2tax reuses config["malt_accession2taxid"] (NCBI nucl_gb.accession2taxid);
##     the aligner's reference names are the FASTA accessions, which it keys on.

if NGSLCA_PROFILERS:

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
        """Assign a lowest-common-ancestor taxid to every read with ngsLCA."""
        output:
            lca="results/NGSLCA/{sample}.lca",
        input:
            bam="results/ALIGNMENT/{sample}.trimmed.namesorted.bam",
            names="results/NGSLCA_DB/names.dmp",
            nodes="results/NGSLCA_DB/nodes.dmp",
            acc2tax=config["malt_accession2taxid"],
        params:
            prefix="results/NGSLCA/{sample}",
            args=config.get("ngslca_args", "-simscorelow 0.85 -simscorehigh 1.0"),
        threads: 4
        log:
            "logs/NGSLCA/{sample}.log",
        conda:
            "../envs/ngslca.yaml"
        benchmark:
            "benchmarks/NGSLCA/{sample}.benchmark.txt"
        message:
            "NgsLCA: ASSIGNING TAXONOMY WITH ngsLCA FOR SAMPLE {input.bam}"
        shell:
            "ngsLCA {params.args} -names {input.names} -nodes {input.nodes} "
            "-acc2tax {input.acc2tax} -bam {input.bam} -outnames {params.prefix} 2> {log}"

    rule NgsLCA_AbundanceMatrix:
        """Count reads per LCA taxid across samples into a wide abundance matrix."""
        output:
            out_dir=directory("results/NGSLCA_ABUNDANCE_MATRIX"),
            abundance_matrix="results/NGSLCA_ABUNDANCE_MATRIX/ngslca_abundance_matrix.txt",
        input:
            lca=expand("results/NGSLCA/{sample}.lca", sample=SAMPLES),
        params:
            exe=WORKFLOW_DIR / "scripts/ngslca_abundance_matrix.py",
            pairs=lambda wildcards, input: " ".join(
                f"{s} {f}" for s, f in zip(SAMPLES, input.lca)
            ),
        threads: 1
        log:
            "logs/NGSLCA_ABUNDANCE_MATRIX/NGSLCA_ABUNDANCE_MATRIX.log",
        conda:
            "../envs/ngslca.yaml"
        message:
            "NgsLCA_AbundanceMatrix: BUILDING ngsLCA ABUNDANCE MATRIX"
        shell:
            "mkdir -p {output.out_dir}; "
            "python {params.exe} {output.abundance_matrix} {params.pairs} 2> {log}"
