#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a wide taxon-abundance matrix from ngsLCA ``.lca`` files.

Each ``.lca`` line is tab-separated: the first field is read metadata
(``readID:seq:len:nhits``) and the second field is the LCA assignment
``taxid:name:rank``; further fields trace the path to root. We count reads per
LCA taxid per sample.

Output (tab-separated), one row per taxid observed in any sample::

    taxid   taxon   <sample1>   <sample2> ...

Usage: ngslca_abundance_matrix.py OUT.tsv SAMPLE1 LCA1 [SAMPLE2 LCA2 ...]
"""
import re
import sys

# taxid:name:rank  (name may itself be quoted and contain spaces/colons)
_LCA = re.compile(r"^(\d+):(.*):([^:]*)$")


def parse_lca(path):
    """Return {taxid: (name, count)} for one .lca file."""
    counts = {}
    names = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            m = _LCA.match(fields[1])
            if not m:
                continue
            taxid = m.group(1)
            name = m.group(2).strip().strip('"')
            if taxid == "1":  # root = unassigned; not a detection
                continue
            counts[taxid] = counts.get(taxid, 0) + 1
            names.setdefault(taxid, name)
    return counts, names


def main() -> None:
    args = sys.argv[1:]
    if len(args) < 3 or (len(args) - 1) % 2 != 0:
        sys.exit("usage: ngslca_abundance_matrix.py OUT.tsv SAMPLE1 LCA1 [SAMPLE2 LCA2 ...]")
    out = args[0]
    pairs = list(zip(args[1::2], args[2::2]))

    samples = []
    per_sample = {}
    names = {}
    for sample, lca in pairs:
        samples.append(sample)
        counts, nm = parse_lca(lca)
        per_sample[sample] = counts
        for t, n in nm.items():
            names.setdefault(t, n)

    taxids = sorted(names, key=int)
    with open(out, "w") as fh:
        fh.write("taxid\ttaxon\t" + "\t".join(samples) + "\n")
        for t in taxids:
            row = [t, names[t]] + [str(per_sample[s].get(t, 0)) for s in samples]
            fh.write("\t".join(row) + "\n")
    print(f"Wrote {len(taxids)} taxa x {len(samples)} samples to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
