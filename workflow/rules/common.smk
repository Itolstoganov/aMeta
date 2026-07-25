import os
import sys
import subprocess as sp
from pathlib import Path
from snakemake.utils import validate, logger, min_version
from snakemake.common import __version__ as snakemake_version
from snakemake.exceptions import WorkflowError
import packaging.version as pv
import pandas as pd
import contextlib
from config import WORKFLOW_DIR
from snakemake.io import Wildcards

min_version("6.0")

# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info("Changing directory from {} to {}".format(CWD, path))

    os.chdir(path)
    try:
        yield
    except Exception as e:
        logger.warning(e)
        logger.warning("Exception caught: ".format(sys.exc_info()[0]))
    finally:
        logger.info("Changing directory back to {}".format(CWD))
        os.chdir(CWD)


##### load config and sample sheets #####
configfile: "config/config.yaml"


if pv.parse(snakemake_version) >= pv.parse("8.0.0"):
    if DeploymentMethod.ENV_MODULES in workflow.deployment_settings.deployment_method:
        envmodules = os.getenv("ANCIENT_MICROBIOME_ENVMODULES", "config/envmodules.yaml")
        configfile: envmodules
else:
    if workflow.use_env_modules:
        envmodules = os.getenv("ANCIENT_MICROBIOME_ENVMODULES", "config/envmodules.yaml")
        configfile: envmodules


validate(config, schema="../schemas/config.schema.yaml")

# mapDamage consumes the Bowtie2 BAM, so it cannot run without bowtie2.
if not config["analyses"].get("bowtie2", True) and config["analyses"].get(
    "mapdamage", True
):
    raise WorkflowError(
        "analyses.bowtie2 is disabled but analyses.mapdamage is enabled; "
        "mapDamage requires the Bowtie2 alignment. Set analyses.mapdamage: "
        "false to skip Bowtie2 alignment."
    )


kw = {"sep": "\t" if config["samplesheet"].endswith(".tsv") else ","}
samples = pd.read_csv(config["samplesheet"], **kw).set_index("sample", drop=False)
samples.index.names = ["sample_id"]
if "exclude" in config["samples"]:
    logger.info(
        "Excluding samples {exclude} from analysis".format(
            exclude=",".join(f"'{x}'" for x in config["samples"]["exclude"])
        )
    )
    samples = samples[~samples.index.isin(config["samples"]["exclude"])]
if "include" in config["samples"]:
    logger.info(
        "Restricting analysis to samples {incl}".format(
            incl=",".join(f"'{x}'" for x in config["samples"]["include"])
        )
    )
    samples = samples[samples.index.isin(config["samples"]["include"])]

validate(samples, schema="../schemas/samples.schema.yaml")

##############################
## Store some workflow metadata
##############################
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_workdir__"] = os.getcwd()
config["__worfklow_commit__"] = None
config["__workflow_commit_link__"] = None

try:
    with cd(workflow.basedir, logger):
        commit = sp.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
        commit_short = (
            sp.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
        )
        config["__workflow_commit__"] = commit_short
        config[
            "__workflow_commit_link__"
        ] = f"https://github.com/NBISweden/aMeta/commit/{commit}"
except Exception as e:
    print(e)
    raise


##############################
# Global variables
##############################
#
# Store some config values in all-caps global vars
#
SAMPLES = samples["sample"].tolist()
# Set the adapter list
ADAPTERS = config["adapters"].get("custom", list())
if config["adapters"]["illumina"]:
    ADAPTERS.append("AGATCGGAAGAG")
if config["adapters"]["nextera"]:
    ADAPTERS.append("CTGTCTCTTATA")

##############################
# Taxonomic profilers
##############################
#
# One or more alignment + taxonomy backends may run together; each writes an
# independent authentication tree under results/AUTHENTICATION/<profiler>/. A
# profiler label is either "malt" (MALT alignment + LCA) or an aligner name
# (aligner + ngsLCA). The scalar `taxonomic_profiler` is honoured for back-compat
# when `profilers` is absent.
PROFILERS = config.get("profilers")
if not PROFILERS:
    _tp = config.get("taxonomic_profiler", "malt")
    if _tp == "malt":
        PROFILERS = ["malt"]
    else:  # legacy "aligner_ngslca": label the profiler by the aligner name
        PROFILERS = [config.get("aligner", {}).get("name", "strobealign")]
# aligner + ngsLCA profilers are every non-malt label.
NGSLCA_PROFILERS = [p for p in PROFILERS if p != "malt"]
# Only a single aligner is configured (one `aligner` block), so every ngsLCA
# profiler label must name that aligner. Running several different aligners with
# ngsLCA simultaneously would additionally need per-aligner ALIGNMENT/NGSLCA paths.
_aligner_name = config.get("aligner", {}).get("name", "strobealign")
_bad = [p for p in NGSLCA_PROFILERS if p != _aligner_name]
if _bad:
    raise WorkflowError(
        f"profilers {_bad} are not 'malt' and do not match the configured "
        f"aligner.name '{_aligner_name}'. Every non-malt profiler is an "
        "aligner+ngsLCA backend and must equal aligner.name; only one aligner "
        "can be configured at a time."
    )


##############################
# Wildcard constraints
##############################
#
# Restrict some globally used wildcards for enhanced performance
#
wildcard_constraints:
    sample=f"({'|'.join(samples['sample'].tolist())})",
    profiler=f"({'|'.join(PROFILERS)})",


##############################
# Input collection functions
##############################
def all_input(wildcards):
    d = {
        "multiqc.after": rules.MultiQC.output,
        "mapdamage": mapdamage_input(wildcards),
        "krakenuniq.krona": krona_input(wildcards),
        "malt.abundance": malt_input(wildcards),
        "auth": authentication_input(wildcards),
        "summary": summary_input(wildcards),
    }
    return d

def get_krakenuniq_preload_option():
    preload_mode = config.get("krakenuniq_preload_mode", "preload-size")
    preload_size = config.get("krakenuniq_preload_size", "32G")

    if preload_mode == "preload":
        if "krakenuniq_preload_size" in config:
            print("Warning: krakenuniq_preload_size is ignored when krakenuniq_preload_mode is set to 'preload'.")
        return "--preload"

    if preload_mode in ["preload-size", "preload_size"]:
        return f"--preload-size {preload_size}"

    raise ValueError("krakenuniq_preload_mode must be 'preload-size', 'preload_size' or 'preload'.")

def mapdamage_input(wildcards):
    if not config["analyses"]["mapdamage"]:
        return []
    return expand("results/MAPDAMAGE/{sample}", sample=SAMPLES)


def authentication_input(wildcards):
    if not config["analyses"]["authentication"]:
        return []
    return expand(
        "results/AUTHENTICATION/{profiler}/.{sample}_done",
        profiler=PROFILERS,
        sample=SAMPLES,
    )


def malt_input(wildcards):
    # Collect the abundance matrix of every configured profiler. Each backend has
    # distinct output paths, so several can be produced in one run.
    out = []
    if "malt" in PROFILERS and config["analyses"]["malt"]:
        out += [
            "results/MALT_ABUNDANCE_MATRIX_SAM/malt_abundance_matrix_sam.txt",
            "results/MALT_ABUNDANCE_MATRIX_RMA6/malt_abundance_matrix_rma6.txt",
        ]
    if NGSLCA_PROFILERS:
        out += ["results/NGSLCA_ABUNDANCE_MATRIX/ngslca_abundance_matrix.txt"]
    return out

def summary_input(wildcards):
    return expand(
        "results/overview_heatmap_scores.{profiler}.pdf", profiler=PROFILERS
    )

def krona_input(wildcards):
    if not config["analyses"]["krona"]:
        return []
    return expand("results/KRAKENUNIQ/{sample}/taxonomy.krona.html", sample=SAMPLES)


def multiqc_input(wildcards):
    """Collect all inputs to multiqc"""
    d = {
        "fastqc_before_trimming": expand(
            "results/FASTQC_BEFORE_TRIMMING/{sample}_fastqc.zip", sample=SAMPLES
        ),
        "fastqc_after_trimming": expand(
            "results/FASTQC_AFTER_TRIMMING/{sample}.trimmed_fastqc.zip",
            sample=SAMPLES,
        ),
        "cutadapt": expand(
            "logs/CUTADAPT_ADAPTER_TRIMMING/{sample}.log", sample=SAMPLES
        ),
    }
    if config["analyses"].get("bowtie2", True):
        d["bowtie2"] = expand("logs/BOWTIE2/{sample}.log", sample=SAMPLES)
    return d


def aggregate_maltextract(wildcards):
    """Collect maltextract output directories"""
    checkpoint_output = checkpoints.Create_Sample_TaxID_Directories.get(
        sample=wildcards.sample, profiler=wildcards.profiler
    ).output[0]
    taxid = glob_wildcards(
        os.path.join(os.path.dirname(checkpoint_output), "{taxid,[0-9]+}")
    ).taxid
    return expand(
        "results/AUTHENTICATION/{profiler}/{sample}/{taxid}/MaltExtract_output/log.txt",
        profiler=wildcards.profiler,
        sample=wildcards.sample,
        taxid=taxid,
    )


def _aggregate_utils(fmt, wildcards):
    """Collect common output for all aggregate functions. Returns a tuple
    of lists sample, and taxid"""
    logger.debug(
        f"Running _aggregate_utils for format '{fmt}', wildcards '{dict(wildcards)}'"
    )
    res = []
    checkpoint_output = checkpoints.Create_Sample_TaxID_Directories.get(
        sample=wildcards.sample, profiler=wildcards.profiler
    ).output[0]
    taxid = glob_wildcards(
        os.path.join(os.path.dirname(checkpoint_output), "{taxid,[0-9]+}")
    ).taxid
    profiler = []
    sample = []
    refid = []
    taxid_out = []
    for tid in taxid:
        wc = Wildcards(
            fromdict={
                "profiler": wildcards.profiler,
                "sample": wildcards.sample,
                "taxid": tid,
            }
        )
        _refid = get_ref_id(wc)
        if _refid is not None and _refid != tid:
            refid.append(_refid)
            taxid_out.append(tid)
            sample.append(wildcards.sample)
            profiler.append(wildcards.profiler)
    if len(refid) > 0:
        res = expand(fmt, zip, profiler=profiler, sample=sample, taxid=taxid_out)
    return res


def aggregate_PMD(wildcards):
    fmt = "results/AUTHENTICATION/{profiler}/{sample}/{taxid}/PMD_plot.frag.pdf"
    return _aggregate_utils(fmt, wildcards)

def aggregate_plots(wildcards):
    fmt = "results/AUTHENTICATION/{profiler}/{sample}/{taxid}/authentic_Sample_{sample}.trimmed.rma6_TaxID_{taxid}.pdf"
    return _aggregate_utils(fmt, wildcards)

def aggregate_scores(wildcards):
    fmt = "results/AUTHENTICATION/{profiler}/{sample}/{taxid}/authentication_scores.txt"
    return _aggregate_utils(fmt, wildcards)

def aggregate_post(wildcards):
    fmt = "results/AUTHENTICATION/{profiler}/{sample}/{taxid}/MaltExtract_output/analysis.RData"
    return _aggregate_utils(fmt, wildcards)

def auth_alignment_sam(wildcards):
    """SAM (accession/tax RNAME space) feeding the authentication path: the MALT
    SAM for the 'malt' profiler, otherwise the aligner's canonical SAM."""
    if wildcards.profiler == "malt":
        return f"results/MALT/{wildcards.sample}.trimmed.sam.gz"
    return f"results/ALIGNMENT/{wildcards.sample}.trimmed.sam.gz"


def get_ref_id(wildcards):
    """Return reference id for a given taxonomy id"""
    ref_id = wildcards.taxid
    infile = f"results/AUTHENTICATION/{wildcards.profiler}/{wildcards.sample}/{wildcards.taxid}/MaltExtract_output/default/readDist/{wildcards.sample}.trimmed.rma6_additionalNodeEntries.txt"
    if not os.path.exists(infile):
        logger.debug(f"No such file {infile}; cannot extract refid")
        return None
    with open(infile) as f:
        contents = f.readlines()
        try:
            ref_id = contents[1].split(";")[1][1:]
        except:
            logger.warning(
                f"Failed to extract ref_id from {infile}; returning taxid {wildcards.taxid}"
            )
            pass
    return ref_id


def format_maltextract_output_directory(wildcards):
    """Format MaltExtract output directory name"""
    return f"results/AUTHENTICATION/{wildcards.profiler}/{wildcards.sample}/{wildcards.taxid}/MaltExtract_output/"
