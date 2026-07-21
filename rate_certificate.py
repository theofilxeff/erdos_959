#!/usr/bin/env python3
"""Executable certificate for Lemma ``Rate certificate``.

The script performs two kinds of checks:

1. Every value k_eta(q) is certified by exact integer inequalities.
2. All transcendental evaluations use mpmath's interval context, and the
   resulting intervals are required to lie strictly inside the bounds printed
   in the manuscript.

Dependency: mpmath >= 1.3.0.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import prod
import sys

import mpmath as mp


ETA_NUM = 659
ETA_DEN = 50_000
T = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
S = (
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 47, 71, 79, 97, 101,
    107, 109, 139, 151, 163, 167, 179,
)

# Fifty decimal digits are far more than the displayed certificate requires.
mp.iv.dps = 50
mp.mp.dps = 100


@dataclass(frozen=True)
class LocalContribution:
    prime: int
    inertia_degree: int
    k: int
    nu: object
    alpha: object


def ramification_index(p: int) -> int:
    return 2 if p == 2 or p in T else 1


def exact_floor_k(q: int) -> int:
    """Return k_eta(q), with the manuscript's two integer certificates."""
    eta = mp.mpf(ETA_NUM) / ETA_DEN
    estimate = int(mp.floor(1 / (mp.power(q, eta) - 1)))

    def first_inequality(k: int) -> bool:
        return (
            pow(q, ETA_NUM) * pow(k, ETA_DEN)
            <= pow(k + 1, ETA_DEN)
        )

    def second_inequality(k: int) -> bool:
        return (
            pow(q, ETA_NUM) * pow(k + 1, ETA_DEN)
            > pow(k + 2, ETA_DEN)
        )

    k = max(0, estimate)
    while not first_inequality(k):
        k -= 1
    while not second_inequality(k):
        k += 1

    assert first_inequality(k)
    assert second_inequality(k)
    return k


def local_contribution(p: int, f: int) -> LocalContribution:
    if f not in (1, 2):
        raise ValueError("The allowed inertia degrees are 1 and 2.")
    e = ramification_index(p)
    k = exact_floor_k(p**f)
    nu = mp.iv.mpf(k) * mp.iv.log(p) / e
    alpha = mp.iv.log(k + 1) / (e * f)
    return LocalContribution(p, f, k, nu, alpha)


def require_interval_inside(name: str, value: object, lower: str, upper: str) -> None:
    lo = mp.iv.mpf(lower)
    hi = mp.iv.mpf(upper)
    if not (value.a > lo and value.b < hi):
        raise AssertionError(
            f"{name} interval {value} is not strictly inside ({lower}, {upper})"
        )


def require_lower_bound(name: str, value: object, bound: str) -> None:
    if not value.a > mp.iv.mpf(bound):
        raise AssertionError(f"{name} has lower endpoint {value.a}, not > {bound}")


def main() -> int:
    eta = mp.iv.mpf(ETA_NUM) / ETA_DEN
    lam = mp.iv.sqrt(4 * prod(T))
    kappa = mp.iv.log(
        2 * mp.iv.sqrt(lam) * mp.iv.log(lam) * mp.iv.e / (4 * mp.iv.pi)
    )

    all_two = [local_contribution(p, 2) for p in S]
    nu_c = sum((item.nu for item in all_two), mp.iv.mpf(0))
    alpha_c = sum((item.alpha for item in all_two), mp.iv.mpf(0))
    rate = -mp.iv.log(2) - eta - kappa + alpha_c
    margin = rate - eta * (
        kappa + mp.iv.log(mp.iv.pi / eta) - mp.iv.log(lam) + nu_c
    )

    # The intervals stated in the appendix.
    require_interval_inside(
        "nu_c", nu_c, "573.7430444817", "573.7430444818"
    )
    require_interval_inside(
        "alpha_c", alpha_c, "19.8003353767", "19.8003353768"
    )
    require_interval_inside(
        "rate", rate, "7.5418334325", "7.5418334327"
    )
    require_interval_inside(
        "margin", margin, "0.0046211747", "0.0046211748"
    )

    # The rounded inequalities used in the main proof.
    require_lower_bound("nu_c", nu_c, "573.743")
    require_lower_bound("rate", rate, "7.541")
    require_lower_bound("margin", margin, "0.0046")

    min_nu = None
    min_alpha = None
    min_margin = None

    for p in S:
        one = local_contribution(p, 1)
        two = local_contribution(p, 2)
        delta_nu = one.nu - two.nu
        delta_alpha = one.alpha - two.alpha
        delta_margin = delta_alpha - eta * delta_nu

        require_lower_bound(f"local nu difference at p={p}", delta_nu, "18.5")
        require_lower_bound(
            f"local alpha difference at p={p}", delta_alpha, "0.96"
        )
        require_lower_bound(
            f"local margin difference at p={p}", delta_margin, "0.70"
        )

        entries = (
            ("nu", delta_nu, min_nu),
            ("alpha", delta_alpha, min_alpha),
            ("margin", delta_margin, min_margin),
        )
        for name, value, current in entries:
            record = (value.a, p, value)
            if current is None or value.a < current[0]:
                if name == "nu":
                    min_nu = record
                elif name == "alpha":
                    min_alpha = record
                else:
                    min_margin = record

    print(f"Python: {sys.version.split()[0]}")
    print(f"mpmath: {mp.__version__}")
    print("All exact floor certificates passed for f=1 and f=2.")
    print(f"nu_c interval:    {nu_c}")
    print(f"alpha_c interval: {alpha_c}")
    print(f"rate interval:    {rate}")
    print(f"margin interval:  {margin}")
    print(
        "Smallest local nu difference:    "
        f"p={min_nu[1]}, interval={min_nu[2]}"
    )
    print(
        "Smallest local alpha difference: "
        f"p={min_alpha[1]}, interval={min_alpha[2]}"
    )
    print(
        "Smallest local margin difference: "
        f"p={min_margin[1]}, interval={min_margin[2]}"
    )
    print("RATE CERTIFICATE VERIFIED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
