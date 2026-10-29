"""Command-line interface for pvalue-adjuster."""
import argparse
import csv
import sys

from .core import METHODS, adjust, qvalue_storey


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pvalue-adjuster", description="Adjust p-values for multiple comparisons.")
    ap.add_argument("pvalues", nargs="*", type=float, help="p-values (or use --csv)")
    ap.add_argument("-m", "--method", action="append", choices=sorted(METHODS), help="method(s); default holm and bh")
    ap.add_argument("--csv", help="CSV file containing a p-value column")
    ap.add_argument("-c", "--column", default="p", help="p-value column name in --csv")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--qvalue", action="store_true", help="also compute Storey q-values")
    ap.add_argument("-o", "--output", help="write CSV with adjusted columns")
    args = ap.parse_args(argv)
    methods = args.method or ["holm", "bh"]
    rows, fields = [], ["p"]
    if args.csv:
        with open(args.csv, encoding="utf-8-sig", newline="") as fh:
            r = csv.DictReader(fh)
            rows, fields = list(r), list(r.fieldnames)
        p = [float(x[args.column]) for x in rows]
    else:
        p = args.pvalues
        rows = [{"p": x} for x in p]
    if not p:
        ap.error("no p-values given")
    for m in methods:
        for row, v in zip(rows, adjust(p, m)):
            row["p_" + m] = "%.6g" % v
        fields.append("p_" + m)
    if args.qvalue:
        q, pi0 = qvalue_storey(p)
        for row, v in zip(rows, q):
            row["q_storey"] = "%.6g" % v
        fields.append("q_storey")
        print("estimated pi0 = %.3f" % pi0, file=sys.stderr)
    fh = open(args.output, "w", newline="", encoding="utf-8") if args.output else sys.stdout
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
    if args.output:
        fh.close()
    for m in methods:
        k = sum(1 for row in rows if float(row["p_" + m]) <= args.alpha)
        print("%-10s %d of %d significant at alpha = %g" % (m, k, len(p), args.alpha), file=sys.stderr)
    return 0
