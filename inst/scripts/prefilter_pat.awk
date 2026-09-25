# Keep only the .pat lines that can overlap a region.
#
# The regions cover a small fraction of the genome, but a .pat.gz is the whole of
# it -- several hundred MB compressed, tens of millions of read patterns.  Reading
# all of that into R to discard 99% of it is the single largest cost in the
# per-sample step, and it is avoidable: both files are sorted by CpG index, so one
# forward pass with a moving pointer is enough.
#
# Usage:
#   zcat sample.pat.gz | awk -v RF=regions.txt -f prefilter_pat.awk > small.pat
#
# regions.txt: two columns, startCpG and endCpG, sorted ascending by startCpG,
# with non-overlapping ranges (what build/index produces).
#
# .pat columns: chr, startCpG, pattern, count.

BEGIN {
    FS = OFS = "\t"
    n = 0
    while ((getline line < RF) > 0) {
        split(line, f, /[ \t]+/)
        if (f[1] ~ /^[0-9]+$/) { n++; rs[n] = f[1] + 0; re[n] = f[2] + 0 }
    }
    close(RF)
    if (n == 0) { print "prefilter: no regions read from " RF > "/dev/stderr"; exit 1 }
    p = 1          # pointer into the region list
    kept = 0; seen = 0
}

{
    seen++
    s = $2 + 0
    e = s + length($3) - 1

    # regions ending before this read starts can never match again
    while (p <= n && re[p] < s) p++
    if (p > n) { nextfile }        # past the last region: nothing further can match

    # the read overlaps region p if it starts at or before that region's end
    if (s <= re[p] && e >= rs[p]) { print; kept++ }
}

END {
    printf("prefilter: kept %d of %d read patterns (%.2f%%)\n",
           kept, seen, seen ? 100.0 * kept / seen : 0) > "/dev/stderr"
}
