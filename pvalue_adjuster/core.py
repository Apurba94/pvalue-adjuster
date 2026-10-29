"""Bonferroni, Sidak, Holm, Hochberg, Hommel-free BH and BY adjustments (same results as R's p.adjust)."""


def _check(p):
    p = [float(x) for x in p]
    if any(not 0 <= x <= 1 for x in p):
        raise ValueError("p-values must be between 0 and 1")
    return p


def bonferroni(p):
    n = len(p)
    return [min(1.0, x * n) for x in p]


def sidak(p):
    n = len(p)
    return [1 - (1 - x) ** n for x in p]


def holm(p):
    n = len(p)
    order = sorted(range(n), key=lambda i: p[i])
    out, running = [0.0] * n, 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (n - rank) * p[i]))
        out[i] = running
    return out


def hochberg(p):
    n = len(p)
    order = sorted(range(n), key=lambda i: p[i], reverse=True)
    out, running = [0.0] * n, 1.0
    for k, i in enumerate(order):          # k = 0 for the largest p
        running = min(running, (k + 1) * p[i])
        out[i] = min(1.0, running)
    return out


def bh(p):
    """Benjamini-Hochberg false discovery rate."""
    n = len(p)
    order = sorted(range(n), key=lambda i: p[i], reverse=True)
    out, running = [0.0] * n, 1.0
    for k, i in enumerate(order):
        rank = n - k
        running = min(running, p[i] * n / rank)
        out[i] = min(1.0, running)
    return out


def by(p):
    """Benjamini-Yekutieli FDR under arbitrary dependence."""
    n = len(p)
    c = sum(1.0 / i for i in range(1, n + 1))
    return [min(1.0, x * c) for x in bh(p)]


METHODS = {"bonferroni": bonferroni, "sidak": sidak, "holm": holm, "hochberg": hochberg, "bh": bh, "fdr": bh,
           "by": by, "none": lambda p: list(p)}


def adjust(pvalues, method="holm"):
    """Return adjusted p-values in the original order. NaN-free floats in [0, 1] required."""
    if method not in METHODS:
        raise ValueError("unknown method %r; choose from %s" % (method, ", ".join(sorted(METHODS))))
    p = _check(pvalues)
    return METHODS[method](p) if p else []


def qvalue_storey(pvalues, lam=0.5):
    """Storey's q-values with pi0 estimated at a single lambda (simple, conservative variant)."""
    p = _check(pvalues)
    n = len(p)
    pi0 = min(1.0, sum(1 for x in p if x > lam) / ((1 - lam) * n)) if n else 1.0
    return [min(1.0, pi0 * q) for q in bh(p)], pi0
