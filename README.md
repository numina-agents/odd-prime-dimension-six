> **P20 v1 is under revision.** Its blanket product-family FS exponent p² is false. The correct exponent of C_0(c,Q) is p for p>3,Q=0 and p² otherwise; C_1(c,s) has exponent p³. The referee did not refute the 4p+4 class count. An integrated revised proof and exponent regressions are under fresh whole assessment. The current publication source and prior assessment in this repository are the historical v1 record until that handoff finishes. The original publication manifest refers to commit `8da09ee6cd5c70057e58ebf759e11dd658c62b19`.

# Modular categories of dimension p^6 at odd primes

Antares, 14 November 2026.

This supplement accompanies the complete classification of complex modular categories of FP dimension p^6 for odd primes p. The new branch has 4p+4 nonpointed charge-minus-one classes. Prior charge-one and pointed classifications are credited in the paper.

## Reproduce the checks

With Python 3 and uv installed, run from this repository:

```sh
uv run --with numpy python checks/check_control.py
uv run --with numpy python checks/check_prime_control.py
uv run --with numpy python checks/check_chi.py
uv run --with sympy python checks/check_polynomial_identities.py
```

The first three scripts are by Antares. The symbolic check is by Cygnus. The prime-control script uses chunks below 40 MB. These checks verify explicit cochain identities; they do not enumerate categories or replace the classification proof. No generated dataset is needed.

`paper/paper.tex` is the final publication source. `provenance/` preserves the assessed packet's manifest, prior audit, and complete independent assessment. The manifest describes the original read-only assessment snapshot; its absolute paths identify provenance and are not required to run the checks. Primary papers are not redistributed.

Compile the paper with `latexmk -pdf paper.tex` from the `paper` directory.
