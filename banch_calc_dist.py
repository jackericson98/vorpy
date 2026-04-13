from __future__ import annotations

import time
import numpy as np

from vorpy.src.calculations.calcs import calc_dist



def calc_dist_py(l0, l1):
    return float(np.sqrt(np.sum(np.square(np.array(l0) - np.array(l1)))))



def time_fn(fn, a, b, n_iter: int) -> float:
    t0 = time.perf_counter()

    for _ in range(n_iter):
        fn(a, b)

    t1 = time.perf_counter()

    return t1 - t0



def main() -> None:
    dim = 3
    n_iter = 2_000_000

    a = np.random.random(dim).astype(np.float64)
    b = np.random.random(dim).astype(np.float64)

    # Warm-up
    calc_dist(a, b)
    calc_dist_py(a, b)

    # Sanity
    print("abs diff:", abs(calc_dist(a, b) - calc_dist_py(a, b)))

    t_accel = time_fn(calc_dist, a, b, n_iter=n_iter)
    t_py = time_fn(calc_dist_py, a, b, n_iter=n_iter)

    print(f"accelerated: {t_accel:.4f}s  ({t_accel / n_iter:.2e}s/call)")
    print(f"python     : {t_py:.4f}s  ({t_py / n_iter:.2e}s/call)")
    print(f"speedup    : {t_py / t_accel:.2f}x")


if __name__ == "__main__":


    main()
