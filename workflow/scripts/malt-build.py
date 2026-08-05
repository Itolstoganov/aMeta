#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2022, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import re
from snakemake.shell import shell
import subprocess as sp

out = sp.run(["malt-build", "--help"], capture_output=True)
regex =  re.compile(r'version (?P<major>\d+)\.(?P<minor>\d+)')
m = regex.search(out.stderr.decode())
if m is None:
    # Assume minor version 4
    minor = 4
else:
    minor = int(m.groupdict()["minor"])

a2t_option = "-a2taxonomy" if minor <= 4 else "-a2t"

log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

accession2taxid = snakemake.params.accession2taxid
max_heap = snakemake.params.max_heap

# The species-subset FASTA is produced by the shared Reference_Subset rule; here we
# only build the MALT database from it.
shell(
    "unset DISPLAY; "
    "malt-build {max_heap} -i {snakemake.input.project_fasta} {a2t_option} {accession2taxid} -s DNA -t {snakemake.threads} -d {snakemake.output.db} {log}"
)
