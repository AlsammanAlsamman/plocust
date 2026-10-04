"""Linkage disequilibrium from genotype matrices."""

from __future__ import annotations

import warnings

import numpy as np


def _standardize(g: np.ndarray) -> np.ndarray:
    """Mean-impute missing calls and scale each column to mean 0, sd 1 (constant columns -> 0)."""
    g = np.array(g, dtype=np.float64, copy=True)
    if g.ndim == 1:
        g = g[:, None]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(g, axis=0)
    mean = np.where(np.isnan(mean), 0.0, mean)
    idx = np.where(np.isnan(g))
    g[idx] = np.take(mean, idx[1])
    g -= g.mean(axis=0)
    sd = g.std(axis=0)
    sd[sd == 0] = np.inf
    return g / sd


def r2_with(g: np.ndarray, target: np.ndarray) -> np.ndarray:
    """r^2 between one variant (`target`, length n) and each column of `g` (n x m)."""
    a = _standardize(target)[:, 0]
    b = _standardize(g)
    return (a @ b / len(a)) ** 2


def r2_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """r^2 between every column of `a` and every column of `b`."""
    sa, sb = _standardize(a), _standardize(b)
    return (sa.T @ sb / sa.shape[0]) ** 2
