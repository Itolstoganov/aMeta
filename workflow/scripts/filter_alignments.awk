# Drop short / heavily soft-clipped local alignments from a SAM stream.
#
# Usage: awk -v min_bp=20 -v min_frac=0.8 -f filter_alignments.awk in.sam

BEGIN { FS = OFS = "\t" }

/^@/ { print; next }

{
    total++
    if (int($2 / 4) % 2 == 1) { unmapped++; print; next }  # 0x4 unmapped

    aligned = 0
    qlen = 0
    cigar = $6
    while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
        n = substr(cigar, RSTART, RLENGTH - 1) + 0
        op = substr(cigar, RSTART + RLENGTH - 1, 1)
        if (op == "M" || op == "I" || op == "=" || op == "X") { aligned += n; qlen += n }
        else if (op == "S" || op == "H") { qlen += n }
        cigar = substr(cigar, RSTART + RLENGTH)
    }

    if (aligned >= min_bp && qlen > 0 && aligned / qlen >= min_frac) { kept++; print }
    else dropped++
}

END {
    printf "filter_alignments: min_bp=%d min_frac=%g records=%d unmapped=%d kept=%d dropped=%d (%.2f%%)\n", \
        min_bp, min_frac, total, unmapped, kept, dropped, \
        (total - unmapped > 0 ? 100 * dropped / (total - unmapped) : 0) > "/dev/stderr"
}
