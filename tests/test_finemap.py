import numpy as np

from plocust.finemap import ld_consistency, susie_rss, susie_rss_multi


def _panel(n, p, rng, twins=()):
    """Genotypes with mild LD; each (a, b) in twins makes column b a copy of column a."""
    base = rng.normal(size=(n, p))
    x = base + 0.6 * np.roll(base, 1, axis=1)  # neighbour correlation
    for a, b in twins:
        x[:, b] = x[:, a]
    return (x - x.mean(0)) / x.std(0)


def _z(x, causal, effects, rng):
    y = x[:, causal] @ np.asarray(effects) + rng.normal(size=len(x))
    y = (y - y.mean()) / y.std()
    return x.T @ y / np.sqrt(len(x))


def test_two_signals_recovered():
    rng = np.random.default_rng(0)
    x = _panel(800, 60, rng)
    z = _z(x, [10, 40], [0.3, 0.25], rng)
    R = np.corrcoef(x, rowvar=False)
    fm = susie_rss(z, R, n=800, L=5)
    assert fm.converged and len(fm.signals) == 2
    found = {int(c) for s in fm.signals for c in s.variants}
    assert {10, 40} <= found
    assert all(len(s.variants) <= 4 for s in fm.signals)


def test_no_signal_gives_no_credible_set():
    rng = np.random.default_rng(1)
    x = _panel(500, 40, rng)
    z = x.T @ rng.normal(size=500) / np.sqrt(500)
    fm = susie_rss(z, np.corrcoef(x, rowvar=False), n=500, L=5)
    assert fm.signals == []


def test_second_population_breaks_ld_ties():
    """Causal SNP 20 has a perfect twin (21) in population 1 only: together the populations pick SNP 20."""
    rng = np.random.default_rng(2)
    x1 = _panel(700, 40, rng, twins=[(20, 21)])
    x2 = _panel(700, 40, rng)
    z1, z2 = _z(x1, [20], [0.3], rng), _z(x2, [20], [0.3], rng)
    R1, R2 = np.corrcoef(x1, rowvar=False), np.corrcoef(x2, rowvar=False)
    single = susie_rss(z1, R1, n=700, L=3)
    assert {20, 21} <= {int(v) for v in single.signals[0].variants}  # cannot tell the twins apart
    both = susie_rss_multi([z1, z2], [R1, R2], [700, 700], L=3)
    top = both.signals[0]
    assert int(top.variants[0]) == 20 and top.pip[0] > 0.9


def test_ld_consistency_flags_flipped_allele():
    rng = np.random.default_rng(3)
    x = _panel(800, 30, rng)
    z = _z(x, [12], [0.35], rng)
    R = np.corrcoef(x, rowvar=False)
    bad = z.copy()
    bad[13] = -bad[13]  # allele flip on a SNP in LD with the causal one
    assert np.argmax(np.abs(ld_consistency(bad, R))) in (12, 13)
    assert np.abs(ld_consistency(bad, R)).max() > np.abs(ld_consistency(z, R)).max()
