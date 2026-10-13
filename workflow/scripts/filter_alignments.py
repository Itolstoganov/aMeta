#!/usr/bin/env python3
"""Drop short or heavily soft-clipped local alignments from a SAM stream.

Usage: filter_alignments.py MIN_ALIGNED_BP MIN_ALIGNED_FRAC < in.sam > out.sam

An alignment is kept when it covers at least MIN_ALIGNED_BP read bases and at
least MIN_ALIGNED_FRAC of the read. Unmapped records pass through.
"""
import sys

import pysam


def main():
    min_bp, min_frac = int(sys.argv[1]), float(sys.argv[2])
    kept = dropped = 0
    with pysam.AlignmentFile("-", "r") as sam_in, pysam.AlignmentFile(
        "-", "w", template=sam_in
    ) as sam_out:
        for rec in sam_in:
            if rec.is_unmapped:
                sam_out.write(rec)
                continue
            # aligned read bases (M/I/=/X) and full read length including clips
            aligned = rec.query_alignment_length
            read_len = rec.infer_read_length() or 0
            if read_len and aligned >= min_bp and aligned >= min_frac * read_len:
                sam_out.write(rec)
                kept += 1
            else:
                dropped += 1
    print(f"filter_alignments: kept {kept}, dropped {dropped}", file=sys.stderr)


if __name__ == "__main__":
    main()
