# Cap the LoF/PC sweep at n-1 PCs

## Goal

Prevent the enrichment pipeline from fitting a residualization model that includes every available PC. If the PC matrix contains `n` PCs, the largest permitted PC count is `max(n - 1, 0)`.

## Behavior

- Adaptive PC grids use `max(n - 1, 0)` as their upper bound.
- Explicit PC-count requests are validated against the same upper bound and reject `n`.
- The reported number of available PCs remains `n`; only the model-selection/residualization limit changes.
- A matrix with zero available PCs retains the intercept-only count `0`; a matrix with one available PC also permits only count `0`.

## Verification

Add regression tests for adaptive and explicit boundaries, including zero/one-PC edge cases, and update existing tests that intentionally used the former full-PC setting.
