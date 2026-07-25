#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Produce MaltExtract-compatible authentication tables from strobealign + ngsLCA.

This is the ngsLCA-path replacement for MEGAN's ``MaltExtract`` (which needs an
``.rma6``). For one (sample, taxid) it reads the aligner's name-sorted BAM plus the
ngsLCA ``.lca`` assignment and writes exactly the files that ``score.R`` and
``authentic.R`` read out of ``MaltExtract_output/`` — but only the first-row values
they actually consume. Everything downstream (score.R, authentic.R,
Breadth_Of_Coverage via get_ref_id) is therefore reused byte-for-byte.

Reads "belong" to the taxid if the taxid appears anywhere in that read's ngsLCA
lineage (i.e. the read resolved to this species or finer). Like MALT, the
per-read metrics (edit distance, damage, read length, identity) are aggregated over
all such member reads at the node level, each contributing its lowest-NM
alignment regardless of which of the species' accessions it landed on. A single
"top reference" (ref_id) -- the target genome the most member reads aligned -- 
is used only for the per-reference outputs (coverage/IGV and the
TOPREFPERCREADS share), matching MaltExtract.

Files written under <out_dir> (= .../MaltExtract_output/), with <B> = rma6 basename:
  log.txt
  default/readDist/<B>_additionalNodeEntries.txt   (encodes ref_id for get_ref_id)
  default/readDist/<B>_alignmentDist.txt            (TotalAlignmentsOnReference, ...)
  default/readDist/<B>_readLengthStat.txt           (Mean, StandardDev)
  default/damageMismatch/<B>_damageMismatch.txt     (C>T_1..20, G>A_1..20, considered)
  default/editDistance/<B>_editDistance.txt         (0..10, higher)
  ancient/editDistance/<B>_editDistance.txt         (0..10, higher; damaged reads only)
  default/percentIdentity/<B>_percentIdentity.txt   (80,85,90,95,100)
  default/filterInformation/<B>_filterTable.txt     (turnedOn?)

Usage:
  ngslca_extract.py --bam BAM --lca LCA --taxid TAXID --node-list NODE_LIST
                    --ref-fasta FASTA --out-dir DIR --rma6-basename B
"""
import argparse
import math
import os
import sys

import pysam

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp_base(b):
    return b.translate(_COMP)


def lca_members(lca_path, taxid):
    """Read IDs whose ngsLCA lineage includes `taxid`.

    The first .lca field is `readID:seq:len:nhits`; readID may itself contain ':'
    (Illumina names), so strip the last three colon-tokens. Remaining tab fields are
    `taxid:name:rank` from the LCA up to root.
    """
    members = set()
    with open(lca_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            parts = fields[0].split(":")
            readid = ":".join(parts[:-3]) if len(parts) > 3 else parts[0]
            for f in fields[1:]:
                if f.split(":", 1)[0] == taxid:
                    members.add(readid)
                    break
    return members


def damage_events(aln, ref_seq):
    """Per matched position: (dist5, dist3, ref_base, read_base) in original-read
    orientation. ref_seq is the forward-strand reference over the aligned span."""
    seq = aln.query_sequence
    if seq is None:
        return []
    qstart = aln.query_alignment_start
    qend = aln.query_alignment_end  # exclusive
    rstart = aln.reference_start
    rev = aln.is_reverse
    out = []
    for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
        rb = ref_seq[rpos - rstart].upper()
        qb = seq[qpos].upper()
        if rev:
            d5 = qend - 1 - qpos
            d3 = qpos - qstart
            rb = revcomp_base(rb)
            qb = revcomp_base(qb)
        else:
            d5 = qpos - qstart
            d3 = qend - 1 - qpos
        out.append((d5, d3, rb, qb))
    return out


class ReadRecord:
    __slots__ = ("ref", "nm", "length", "aln_len", "events")

    def __init__(self, ref, nm, length, aln_len, events):
        self.ref = ref
        self.nm = nm
        self.length = length
        self.aln_len = aln_len
        self.events = events


def scan_bam(bam_path, ref_fasta, members):
    """Scan the member reads' alignments.

    Returns (best, any_on_ref):
      best       = {readid: ReadRecord} of each member's best (lowest-NM) alignment.
      any_on_ref = {ref: set(readid)} of every member read that has any alignment on
                   a reference, used only to pick the top reference and its
                   TOPREFPERCREADS share.
    """
    best = {}
    any_on_ref = {}
    fasta = pysam.FastaFile(ref_fasta)
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for aln in bam:
            if aln.is_unmapped or aln.query_name not in members:
                continue
            qname = aln.query_name
            ref = aln.reference_name
            any_on_ref.setdefault(ref, set()).add(qname)
            try:
                nm = aln.get_tag("NM")
            except KeyError:
                nm = None
            ref_seq = fasta.fetch(ref, aln.reference_start, aln.reference_end)
            events = damage_events(aln, ref_seq)
            if nm is None:
                nm = sum(1 for _, _, rb, qb in events if rb != qb)
            length = aln.infer_read_length() or aln.query_length or len(events)
            aln_len = aln.query_alignment_length or len(events)
            prev = best.get(qname)
            if prev is None or nm < prev.nm:
                best[qname] = ReadRecord(ref, nm, length, aln_len, events)
    fasta.close()
    return best, any_on_ref


def is_ancient(events):
    """Terminal-damage read: C>T at the 5' terminus or G>A at the 3' terminus."""
    for d5, d3, rb, qb in events:
        if d5 == 0 and rb == "C" and qb == "T":
            return True
        if d3 == 0 and rb == "G" and qb == "A":
            return True
    return False


def write_table(path, header_cols, rowname, values):
    """Whitespace/tab table: header (no rowname col) then one data row."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(header_cols) + "\n")
        fh.write(rowname + "\t" + "\t".join(str(v) for v in values) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--lca", required=True)
    ap.add_argument("--taxid", required=True)
    ap.add_argument("--node-list", required=True)
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--rma6-basename", required=True)
    args = ap.parse_args()

    taxid = args.taxid
    B = args.rma6_basename
    out = args.out_dir
    rd_dir = os.path.join(out, "default", "readDist")
    log_path = os.path.join(out, "log.txt")

    members = lca_members(args.lca, taxid)
    best, any_on_ref = scan_bam(args.bam, args.ref_fasta, members) if members else ({}, {})

    # Node-level read set: every member read that aligned, one record (its best
    # alignment) each. 
    node_recs = list(best.values())
    n = len(node_recs)

    best_refs = {rec.ref for rec in node_recs}
    if best_refs:
        # most member reads aligned anywhere; ref name breaks ties deterministically
        ref_id = max(best_refs, key=lambda r: (len(any_on_ref.get(r, ())), r))
        top_ref_reads = len(any_on_ref.get(ref_id, ()))
    else:
        ref_id = taxid
        top_ref_reads = 0

    encoded_ref = ref_id if best_refs else taxid
    pct = round(100.0 * top_ref_reads / n, 1) if n else 0.0
    os.makedirs(rd_dir, exist_ok=True)
    with open(os.path.join(rd_dir, f"{B}_additionalNodeEntries.txt"), "w") as fh:
        fh.write("Node\tadditionalNodeEntries\n")
        fh.write(f"{taxid}\t{pct};_{encoded_ref};_TOPREFPERCREADS\n")

    # alignmentDist
    write_table(
        os.path.join(rd_dir, f"{B}_alignmentDist.txt"),
        ["TotalAlignmentsOnReference", "nonDuplicatesonReference",
         "uniquePerReference", "nonStacked"],
        taxid, [n, n, n, n],
    )

    # readLengthStat
    lengths = [rec.length for rec in node_recs]
    mean = sum(lengths) / n if n else 0
    sd = math.sqrt(sum((x - mean) ** 2 for x in lengths) / n) if n else 0
    write_table(
        os.path.join(rd_dir, f"{B}_readLengthStat.txt"),
        ["Mean", "StandardDev"], taxid, [round(mean, 4), round(sd, 4)],
    )

    # damageMismatch: C>T_1..20 (5' 1-10, 3' -10..-1), G>A_1..20, then considered
    ct5 = [0] * 10; c5 = [0] * 10; ct3 = [0] * 10; c3 = [0] * 10
    ga5 = [0] * 10; g5 = [0] * 10; ga3 = [0] * 10; g3 = [0] * 10
    for rec in node_recs:
        for d5, d3, rb, qb in rec.events:
            if d5 < 10:
                if rb == "C":
                    c5[d5] += 1
                    if qb == "T":
                        ct5[d5] += 1
                if rb == "G":
                    g5[d5] += 1
                    if qb == "A":
                        ga5[d5] += 1
            if d3 < 10:
                if rb == "C":
                    c3[d3] += 1
                    if qb == "T":
                        ct3[d3] += 1
                if rb == "G":
                    g3[d3] += 1
                    if qb == "A":
                        ga3[d3] += 1

    ct = []
    for i in range(10):                       # C>T_1..C>T_10  (5')
        ct.append(round(ct5[i] / c5[i], 6) if c5[i] else 0.0)
    for i in range(9, -1, -1):                # C>T_11..C>T_20 (3', 20 = -1)
        ct.append(round(ct3[i] / c3[i], 6) if c3[i] else 0.0)
    for i in range(10):                       # G>A_1..G>A_10  (5')
        ct.append(round(ga5[i] / g5[i], 6) if g5[i] else 0.0)
    for i in range(9, -1, -1):                # G>A_11..G>A_20 (3', 20 = -1)
        ct.append(round(ga3[i] / g3[i], 6) if g3[i] else 0.0)
    dam_cols = ([f"C>T_{i}" for i in range(1, 21)]
                + [f"G>A_{i}" for i in range(1, 21)] + ["considered"])
    write_table(
        os.path.join(out, "default", "damageMismatch", f"{B}_damageMismatch.txt"),
        dam_cols, taxid, ct + [n],
    )

    # editDistance (all + ancient): bins 0..10, higher
    def edit_row(records):
        bins = [0] * 11
        higher = 0
        for rec in records:
            if rec.nm <= 10:
                bins[rec.nm] += 1
            else:
                higher += 1
        return bins + [higher]

    edit_cols = [str(i) for i in range(11)] + ["higher"]
    write_table(
        os.path.join(out, "default", "editDistance", f"{B}_editDistance.txt"),
        edit_cols, taxid, edit_row(node_recs),
    )
    ancient = [rec for rec in node_recs if is_ancient(rec.events)]
    write_table(
        os.path.join(out, "ancient", "editDistance", f"{B}_editDistance.txt"),
        edit_cols, taxid, edit_row(ancient),
    )

    # percentIdentity: bins 80,85,90,95,100
    pid_bins = {"80": 0, "85": 0, "90": 0, "95": 0, "100": 0}
    for rec in node_recs:
        ident = 100.0 * (1 - rec.nm / rec.aln_len) if rec.aln_len else 0.0
        if ident >= 97.5:
            pid_bins["100"] += 1
        elif ident >= 92.5:
            pid_bins["95"] += 1
        elif ident >= 87.5:
            pid_bins["90"] += 1
        elif ident >= 82.5:
            pid_bins["85"] += 1
        else:
            pid_bins["80"] += 1
    write_table(
        os.path.join(out, "default", "percentIdentity", f"{B}_percentIdentity.txt"),
        ["80", "85", "90", "95", "100"], taxid,
        [pid_bins[k] for k in ("80", "85", "90", "95", "100")],
    )

    # filterTable
    write_table(
        os.path.join(out, "default", "filterInformation", f"{B}_filterTable.txt"),
        ["turnedOn?"], taxid, [0],
    )

    with open(log_path, "w") as fh:
        fh.write(
            f"ngslca_extract: taxid={taxid} members={len(members)} node_reads={n} "
            f"ref_id={ref_id} top_ref_reads={top_ref_reads} ({pct}%) "
            f"ancient={len(ancient) if n else 0}\n"
        )
    print(f"taxid {taxid}: {len(members)} members, {n} node reads, "
          f"ref_id {ref_id} ({pct}% of node)", file=sys.stderr)


if __name__ == "__main__":
    main()
