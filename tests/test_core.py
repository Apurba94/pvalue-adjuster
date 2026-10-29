import io
import os
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

from pvalue_adjuster.cli import main
from pvalue_adjuster.core import adjust, qvalue_storey

P = [0.01, 0.02, 0.03, 0.04, 0.05]
# Values checked against R: p.adjust(c(0.01, 0.04, 0.03, 0.005, 0.2), method = ...)
Q = [0.01, 0.04, 0.03, 0.005, 0.2]


def r(xs, n=6):
    return [round(x, n) for x in xs]


class AdjustTests(unittest.TestCase):
    def test_holm_and_hochberg(self):
        self.assertEqual(r(adjust(P, "holm")), [0.05, 0.08, 0.09, 0.09, 0.09])
        self.assertEqual(r(adjust(P, "hochberg")), [0.05] * 5)

    def test_bh_and_by(self):
        self.assertEqual(r(adjust(P, "bh")), [0.05] * 5)
        self.assertEqual(r(adjust(P, "by")), [round(0.05 * (1 + 1 / 2 + 1 / 3 + 1 / 4 + 1 / 5), 6)] * 5)

    def test_against_r_unsorted(self):
        self.assertEqual(r(adjust(Q, "bonferroni")), [0.05, 0.2, 0.15, 0.025, 1.0])
        self.assertEqual(r(adjust(Q, "holm")), [0.04, 0.09, 0.09, 0.025, 0.2])
        self.assertEqual(r(adjust(Q, "bh")), [0.025, 0.05, 0.05, 0.025, 0.2])

    def test_sidak_and_bounds(self):
        self.assertAlmostEqual(adjust([0.01, 0.5], "sidak")[0], 1 - 0.99 ** 2)
        with self.assertRaises(ValueError):
            adjust([1.2])
        with self.assertRaises(ValueError):
            adjust([0.1], "magic")
        self.assertEqual(adjust([], "bh"), [])

    def test_qvalue(self):
        p = [0.001, 0.002, 0.003, 0.004, 0.01, 0.02, 0.6, 0.9]
        q, pi0 = qvalue_storey(p)
        self.assertAlmostEqual(pi0, 2 / (0.5 * 8))          # two p-values above lambda = 0.5
        self.assertEqual([round(x, 6) for x in q], [round(0.5 * x, 6) for x in adjust(p, "bh")])
        self.assertEqual(qvalue_storey([0.6, 0.7, 0.8])[1], 1.0)   # pi0 is capped at 1

    def test_cli_csv(self):
        with tempfile.TemporaryDirectory() as d:
            src, out = os.path.join(d, "p.csv"), os.path.join(d, "adj.csv")
            with open(src, "w", encoding="utf-8") as fh:
                fh.write("gene,pval\n" + "\n".join("g%d,%s" % (i, p) for i, p in enumerate(Q)) + "\n")
            with redirect_stderr(io.StringIO()) as err, redirect_stdout(io.StringIO()):
                main(["--csv", src, "-c", "pval", "-m", "bh", "-o", out])
            self.assertIn("bh         4 of 5 significant", err.getvalue())
            with open(out, encoding="utf-8") as fh:
                self.assertIn("g0,0.01,0.025", fh.read())


if __name__ == "__main__":
    unittest.main()
