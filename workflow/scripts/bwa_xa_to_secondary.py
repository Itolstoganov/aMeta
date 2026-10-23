#!/usr/bin/env python3
"""Turn the XA alternative hits of bwa samse into secondary SAM records.

Usage: bwa_xa_to_secondary.py MAX_HITS < in.sam > out.sam

ngsLCA needs every hit of a read as its own record. bwa samse lists them only in
the XA tag (chr,+-pos,CIGAR,NM), and omits XA altogether when a read has more
than samse -n (MAX_HITS) best and suboptimal hits (X0 + X1). Such reads are
dropped: their primary alone would give a falsely specific LCA. XA hits running
past the end of their reference, which bwa reports for reads bridging two
concatenated references, are dropped as well. Secondary records carry SEQ and
QUAL, reverse-complemented when the hit is on the other strand, so that damage
can be read from them.
"""
import re
import sys

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def ref_span(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cigar))


def main():
    max_hits = int(sys.argv[1])
    primaries = secondaries = dropped = off_end = 0
    with pysam.AlignmentFile("-", "r") as sam_in, pysam.AlignmentFile(
        "-", "w", template=sam_in
    ) as sam_out:
        for rec in sam_in:
            if rec.is_unmapped:
                sam_out.write(rec)
                continue
            n_hits = (rec.get_tag("X0") if rec.has_tag("X0") else 1) + (
                rec.get_tag("X1") if rec.has_tag("X1") else 0
            )
            if not rec.has_tag("XA") and n_hits > max_hits:
                dropped += 1
                continue
            alts = rec.get_tag("XA").rstrip(";").split(";") if rec.has_tag("XA") else []
            rec.set_tag("XA", None)
            sam_out.write(rec)
            primaries += 1
            seq, qual = rec.query_sequence, rec.query_qualities
            for alt in alts:
                chrom, pos, cigar, nm = alt.split(",")
                start = abs(int(pos)) - 1
                if start + ref_span(cigar) > sam_out.header.get_reference_length(chrom):
                    off_end += 1
                    continue
                sec = pysam.AlignedSegment(sam_out.header)
                sec.query_name = rec.query_name
                sec.flag = 0x100 | (0x10 if pos[0] == "-" else 0)
                sec.reference_name = chrom
                sec.reference_start = start
                sec.mapping_quality = 0
                sec.cigarstring = cigar
                if sec.is_reverse == rec.is_reverse:
                    sec.query_sequence = seq
                    sec.query_qualities = qual
                else:
                    sec.query_sequence = seq.translate(COMPLEMENT)[::-1]
                    sec.query_qualities = qual[::-1] if qual is not None else None
                sec.set_tag("NM", int(nm))
                sam_out.write(sec)
                secondaries += 1
    print(
        f"bwa_xa_to_secondary: {primaries} primary, {secondaries} secondary records; "
        f"dropped {dropped} reads with more than {max_hits} hits, "
        f"{off_end} hits running past the reference end",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
