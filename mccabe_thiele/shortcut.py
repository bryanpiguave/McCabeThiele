"""Fenske-Underwood-Gilliland shortcut method (binary), used to cross-check the
graphical McCabe-Thiele result."""

import math
from scipy.optimize import brentq


def fenske(x_D, x_W, alpha):
    """Minimum number of equilibrium stages (including the reboiler) at total reflux."""
    return math.log((x_D / (1 - x_D)) * ((1 - x_W) / x_W)) / math.log(alpha)


def underwood(x_D, z, q, alpha):
    """Minimum reflux ratio for a binary with constant relative volatility."""
    def g(theta):
        return alpha * z / (alpha - theta) + (1 - z) / (1 - theta) - (1 - q)
    theta = brentq(g, 1 + 1e-9, alpha - 1e-9)
    return alpha * x_D / (alpha - theta) + (1 - x_D) / (1 - theta) - 1


def gilliland(R, R_min, N_min):
    """Number of equilibrium stages from the Molokanov form of the Gilliland correlation."""
    X = (R - R_min) / (R + 1)
    Y = 1 - math.exp((1 + 54.4 * X) / (11 + 117.2 * X) * (X - 1) / math.sqrt(X))
    return (Y + N_min) / (1 - Y)


def kirkbride(N, x_D, x_W, z, D, W):
    """Feed stage (from the top) from the Kirkbride equation."""
    ratio = ((1 - z) / z * (x_W / (1 - x_D)) ** 2 * W / D) ** 0.206   # N_rect / N_strip
    return N * ratio / (1 + ratio)


def fug(vle, x_D, x_W, feed, R, D, W):
    alpha = math.sqrt(float(vle.alpha(x_D)) * float(vle.alpha(x_W)))
    N_min = fenske(x_D, x_W, alpha)
    R_min = underwood(x_D, feed.z, feed.q, alpha)
    N = gilliland(R, R_min, N_min) if R > R_min else float("inf")
    return {"alpha": alpha, "N_min": N_min, "R_min": R_min, "N": N,
            "feed_stage": kirkbride(N, x_D, x_W, feed.z, D, W)}
