# [Catalan's constant is irrational](paper.pdf)

**Author:** OpenAI

**Date:** September 24, 2026

## Citation

```bibtex
@misc{OAI:Catalans-constant-is-irrational-September-24-2026,
  author = {{OpenAI}},
  title = {{Catalan's constant is irrational}},
  howpublished = {OpenAI Math Release preprint
                  \href{https://github.com/openai/math/blob/main/preprints/Catalans-constant-is-irrational-September-24-2026/paper.pdf}{OAI:Catalans-constant-is-irrational-September-24-2026}},
  year = {2026}
}
```

## Verify the fixed-matrix certificate

Run `python3 verify_fixed_matrices.py` from this directory. The script uses only the Python standard library and reconstructs the three 48-by-48 matrices modulo 101 from Section 4.4. It checks all 144 pivots and the stated row swaps against Table 1, and raises an error if they differ.

This verifies the fixed-matrix arithmetic certificate only. It does not verify the growing-matrix reduction, the analytic estimates, or the full irrationality theorem. Lean verification is described in the repository’s `lean/` documentation.
