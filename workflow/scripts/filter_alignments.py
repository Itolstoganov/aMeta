#!/usr/bin/env python3
"""Drop short or heavily soft-clipped local alignments from a SAM stream.

Usage: filter_alignments.py MIN_ALIGNED_BP MIN_ALIGNED_FRAC [MAX_SCORE_DROP] < in.sam > out.sam

An alignment is kept when it covers at least MIN_ALIGNED_BP read bases and at
least MIN_ALIGNED_FRAC of the read. With MAX_SCORE_DROP, an alignment is also
dropped when its AS tag is more than MAX_SCORE_DROP below the best AS of the
read; this needs the read's records to be consecutive, as aligners emit them.
Unmapped records pass through.
"""
import itertools
import sys

import pysam


def main():
    min_bp, min_frac = int(sys.argv[1]), float(sys.argv[2])
    max_drop = float(sys.argv[3]) if len(sys.argv) > 3 else None
    kept = dropped = 0
    with pysam.AlignmentFile("-", "r") as sam_in, pysam.AlignmentFile(
        "-", "w", template=sam_in
    ) as sam_out:
        for _, group in itertools.groupby(sam_in, key=lambda rec: rec.query_name):
            group = list(group)
            min_score = None
            if max_drop is not None:
                scores = [rec.get_tag("AS") for rec in group if not rec.is_unmapped]
                if scores:
                    min_score = max(scores) - max_drop
            for rec in group:
                if rec.is_unmapped:
                    sam_out.write(rec)
                    continue
                # aligned read bases (M/I/=/X) and full read length including clips
                aligned = rec.query_alignment_length
                read_len = rec.infer_read_length() or 0
                if (
                    read_len
                    and aligned >= min_bp
                    and aligned >= min_frac * read_len
                    and (min_score is None or rec.get_tag("AS") >= min_score)
                ):
                    sam_out.write(rec)
                    kept += 1
                else:
                    dropped += 1
    print(f"filter_alignments: kept {kept}, dropped {dropped}", file=sys.stderr)


if __name__ == "__main__":
    main()
