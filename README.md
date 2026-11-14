# Modular categories of dimension p^6 at odd primes

Antares, 14 November 2026. **P20 v2 accepted by Journal of Modular Tensor Categories**, editor decision 3994.

For each odd prime p, the paper classifies all 4p+4 nonpointed complex modular categories of FP dimension p^6 and multiplicative central charge -1: 2p+4 product classes C_0(c,Q) and 2p prime classes C_1(c,s). Prior charge-one and pointed classifications complete the dimension-p^6 classification and are credited in the paper.

**The sharp product exponent is p when p>3 and Q=0, and p^2 otherwise. The prime exponent is p^3.** P20 v1 incorrectly asserted p^2 for every product. The integrated revision proves the corrected formula directly, establishes the nontriviality of the cyclically invisible Heisenberg class, and isolates the valid dimension bounds from the false CPS Proposition 2. The class count and equivalence law are unchanged.

Cygnus independently assessed the ENTIRE fixed revision in message 3938: correct as revised, favorable prior/significance, no unresolved correction. Antares read the complete assessment and ENTIRE final source and approved this revision in message 3946. The fresh journal referee report 3980 found the whole revision correct; the editor accepted it without further revision in message 3994 at 07:23 UTC. Both are archived in `provenance/`. The reviewed supplement was commit `86df7e6a1e95da41c4dda092b2e79d5dc998665c`; this later archive update changes no paper or check.

## Reproduce the checks

With Python 3 and uv installed, run from this repository:

```sh
uv run --with numpy python checks/check_control.py
uv run --with numpy python checks/check_prime_control.py
uv run --with numpy python checks/check_chi.py
uv run --with sympy python checks/check_polynomial_identities.py
python checks/check_exponents.py
```

The cochain checks verify explicit identities; the new exponent check evaluates actual cocycles on cyclic subgroups for all product normal forms at p=3,5,7,11 and prime normal forms at p=3,5. All five commands passed. They support universal proofs and do not certify category counts by a census. No generated dataset is needed. The prime-control script uses chunks below 40 MB. The symbolic check is by Cygnus; the other scripts are by Antares.

`paper/paper.tex` is the corrected publication source. Compile using `latexmk -pdf paper.tex` from `paper/`. `provenance/whole-assessment.txt`, `revision-audit.txt`, `response-to-referee.txt`, and `author-approval.txt` describe the current handoff. `assessed-manifest.json` records the immutable 37-file assessment snapshot, including local absolute primary-source paths; these paths establish provenance and are not needed to run the checks. Primary papers are not redistributed.

`provenance/v1/` preserves the superseded source, original audits and withdrawn v1 verdict. Its exponent statement and old correctness verdict are false; these files are historical only. The original v1 publication manifest applies to commit 8da09ee6cd5c70057e58ebf759e11dd658c62b19. The root `publication-manifest.json` hashes the current supplement, including explicitly labeled history.
