"""
Standalone validation for the symmetric-Frobenius projection.

Runs the dual root-finder against the CVXPY reference solver on 20 random
n=10 problems and prints PASS / FAIL with the worst-case KKT residuals.
Exits 0 on PASS, 1 on FAIL. Intended as the first sanity check after a
fresh clone + pip install.

Usage:
    python validate_projection.py
"""

import sys

try:
    import phase2_benchmark as p2
except ImportError as exc:
    print(f"FAIL: could not import phase2_benchmark ({exc}).", file=sys.stderr)
    print("Did you run `pip install -r requirements.txt`?", file=sys.stderr)
    sys.exit(2)


def main() -> int:
    try:
        result = p2.unit_test_projection(n_problems=20, n=10)
    except RuntimeError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"FAIL: unexpected error during projection test: {exc!r}", file=sys.stderr)
        return 1

    print(f"PASS: max ||X_dual - X_cvxpy||_F = {result['max_diff']:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
