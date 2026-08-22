#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2022, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
import struct
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

# two-column, version-less accession -> taxid, prepared by rule Malt_Acc2Taxa
accession2taxid = snakemake.input.acc2taxa
max_heap = snakemake.params.max_heap
# passed explicitly in both directions so the build log always records which
# labelling produced the database (malt-build's own default is true)
parse_taxon_names = "true" if snakemake.params.parse_taxon_names else "false"


def check_a2t_format(path):
    """Reject an -a2t file malt-build would read as all-zero, before building.

    malt-build wants '<version-less accession><TAB><taxid>' and silently stores
    taxon 0 for anything else, so a 41-minute build would otherwise be needed to
    discover the mistake. The post-build check below is the real backstop; this
    only makes the common case fail immediately.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1].isdigit():
                return
            raise RuntimeError(
                "{} is not in malt-build's -a2t format: expected "
                "'<version-less accession><TAB><taxid>', got {!r}. NCBI's "
                "four-column accession2taxid needs converting first (see rule "
                "Malt_Acc2Taxa).".format(path, line.rstrip("\n")[:120])
            )
    raise RuntimeError("{} is empty; malt-build would label no references".format(path))


check_a2t_format(accession2taxid)

# The species-subset FASTA is produced by the shared Reference_Subset rule; here we
# only build the MALT database from it.
shell(
    "unset DISPLAY; "
    "malt-build {max_heap} -i {snakemake.input.project_fasta} {a2t_option} {accession2taxid} -s DNA -t {snakemake.threads} -d {snakemake.output.db} -tn {parse_taxon_names} {log}"
)


def count_mapped_references(db_dir):
    """Count how many references the finished database labels with a taxon.

    taxonomy.idx is a printable magic string ("MATaxonomyV1.1"), a big-endian
    int32 reference count, then one big-endian int32 taxon id per reference. A
    reference malt-build could not label is stored as taxon 0. Returns
    (mapped, total), or (None, None) if the file does not parse -- a future MALT
    may change the layout, and that should not fail an otherwise good build.
    """
    path = os.path.join(db_dir, "taxonomy.idx")
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError:
        return None, None
    if not data.startswith(b"MATaxonomy"):
        return None, None
    # the version suffix on the magic varies, so measure the printable run
    offset = 0
    while offset < len(data) and 32 <= data[offset] < 127:
        offset += 1
    try:
        (total,) = struct.unpack_from(">i", data, offset)
        taxids = struct.unpack_from(">{}i".format(total), data, offset + 4)
    except (struct.error, ValueError):
        return None, None
    return sum(1 for taxid in taxids if taxid != 0), total


def report(message):
    print(message)
    with open(snakemake.log[0], "a") as fh:
        print(message, file=fh)


# malt-build maps every reference to taxon 0 rather than failing when it cannot
# read the -a2t file (a wrong column layout, say). malt-run then aligns normally
# and bins nothing -- "Assig. Taxonomy: 0" -- and the empty classification only
# surfaces much later, as MaltExtract tables full of NA. Catch it here, while the
# cause is still obvious.
mapped, total = count_mapped_references(snakemake.output.db)
if mapped is None:
    report(
        "WARNING: could not read taxonomy.idx; skipped the taxon-mapping check"
    )
elif mapped == 0:
    raise RuntimeError(
        "malt-build labelled 0 of {:,} references with a taxon, so this database "
        "would classify no reads at all. Check that {} is a two-column "
        "'<version-less accession><TAB><taxid>' file{}.".format(
            total,
            accession2taxid,
            "" if parse_taxon_names == "true"
            else ", since malt_parse_taxon_names is false and the accession "
                 "mapping is the only source of labels",
        )
    )
elif mapped < total:
    report(
        "malt-build labelled {:,} of {:,} references with a taxon "
        "({:,} unmapped)".format(mapped, total, total - mapped)
    )
else:
    report("malt-build labelled all {:,} references with a taxon".format(total))
