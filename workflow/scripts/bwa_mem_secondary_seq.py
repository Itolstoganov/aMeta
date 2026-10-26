#!/usr/bin/env python3
"""Give the secondary records of bwa mem -a the SEQ and QUAL of their read.

Usage: bwa_mem_secondary_seq.py < in.sam > out.sam

bwa mem writes SEQ and QUAL as '*' on secondary records, but damage is read from
every record a read is assigned through. They are copied from the read's primary
record, reverse-complemented when the secondary hit is on the other strand and
trimmed by its hard clips, if any. This needs the read's records to be
consecutive, as bwa mem emits them.
"""
import itertools
import sys

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def hard_clips(rec):
    """(left, right) hard-clipped read bases of a record."""
    cigar = rec.cigartuples or []
    left = cigar[0][1] if cigar and cigar[0][0] == pysam.CHARD_CLIP else 0
    right = cigar[-1][1] if len(cigar) > 1 and cigar[-1][0] == pysam.CHARD_CLIP else 0
    return left, right


def main():
    filled = missing = 0
    with pysam.AlignmentFile("-", "r") as sam_in, pysam.AlignmentFile(
        "-", "w", template=sam_in
    ) as sam_out:
        for _, group in itertools.groupby(sam_in, key=lambda rec: rec.query_name):
            group = list(group)
            primary = next(
                (r for r in group if not r.is_secondary and not r.is_supplementary
                 and r.query_sequence), None)
            for rec in group:
                if rec.is_secondary and not rec.query_sequence:
                    if primary is None:
                        missing += 1
                        continue
                    seq, qual = primary.query_sequence, primary.query_qualities
                    if rec.is_reverse != primary.is_reverse:
                        seq = seq.translate(COMPLEMENT)[::-1]
                        qual = qual[::-1] if qual is not None else None
                    left, right = hard_clips(rec)
                    end = len(seq) - right
                    # setting SEQ resets QUAL, so QUAL goes second
                    rec.query_sequence = seq[left:end]
                    rec.query_qualities = qual[left:end] if qual is not None else None
                    filled += 1
                sam_out.write(rec)
    print(
        f"bwa_mem_secondary_seq: filled {filled} secondary records, "
        f"dropped {missing} without a primary record",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
