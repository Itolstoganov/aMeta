rule Plot_Authentication_Score:
    output:
        heatmap="results/overview_heatmap_scores.{profiler}.pdf",
    input:
        scores=expand(
            "results/AUTHENTICATION/{{profiler}}/.{sample}_done", sample=SAMPLES
        ),
    message:
        "Plot_Authentication_Score: PLOTTING HEATMAP OF AUTHENTICATION SCORES FOR PROFILER {wildcards.profiler}"
    params:
        exe=WORKFLOW_DIR / "scripts/plot_score.R",
    log:
        "logs/PLOT_AUTHENTICATION_SCORE/{profiler}.log",
    threads: 1
    conda:
        "../envs/r.yaml"
    envmodules:
        *config["envmodules"]["r"],
    shell:
        "Rscript {params.exe} results/AUTHENTICATION/{wildcards.profiler} results "
        "overview_heatmap_scores.{wildcards.profiler} &> {log}"

