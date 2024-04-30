"""
pyCHomP Benchmarks
==================
Run a series of benchmarks and logs the results in the command line
and/or a file, as well as generating plots.

Global and benchmark-specific options can be configured at the top, and
runs can be added/removed/customized in each relevant section.
"""

import timeit
import os
import matplotlib.pyplot as plt
from pychomp import *


"""
Global Configuration Options
----------------------------
cutoff_time: float
    If a run takes this long, skip to the next method instead of the
    next run. Can be set to `float("Inf")` to avoid this behavior.
iterations: int
    The number of times each benchmark is run for. Results are averaged.
plot_filename: str, optional
    Filename of the plot within the `pyCHomP2/examples/plots` directory.
    If not set, the plot will not be generated.
report_filename: str, optional
    Filename of the report file within the `pyCHomP2/examples/reports`
    directory. If not set, times will not be saved.
verbose: bool
    Whether to print to command line.
progress: bool
    Whether to print progress bars for homology computations
"""
cutoff_time = 20.0
iterations = 20
plot_filename = "benchmark.png"
report_filename = "benchmark.txt"
verbose = True
progress = True


def benchmark(
    name,
    iterations,
    method_configs,
    run_configs,
    stmt_function,
    setup_function=(lambda **kw: None),
    verbose=True,
    file=None
):
    """
    Runs a specified benchmark using the `timeit` package.

    Parameters
    ----------
    name: str
    iterations: int
    method_configs, run_configs: dict
        See examples of config formatting below. These are arguments
        passed to the `stmt_function` and `setup_function` parameters.
    stmt_function, setup_function: function
        Passed to the `timeit` function.
    verbose: bool
        Whether to print to command line.
    file: file, optional
        Logging file.

    Returns
    -------
    dict
        Dictionary with methods as keys and lists of times as values.
    """

    if verbose:
        print("*" * (len(name) + 4))
        print(f"* {name} *")
        print("*" * (len(name) + 4))
    if file is not None:
        print("*" * (len(name) + 4), file=file)
        print(f"* {name} *", file=file)
        print("*" * (len(name) + 4), file=file)

    times_dict = {} # Stores times for return
    method_counter = 1
    num_method = len(method_configs)
    num_run = len(run_configs)

    # Complete all runs for each method before iterating to next method.
    for method, method_config in method_configs.items():
        times = []
        run_counter = 1

        method_string = f"{method} ({method_counter}/{num_method})"
        if verbose:
            print(method_string)
            print("-" * len(method_string))
        if file is not None:
            print(method_string, file=file)
            print("-" * len(method_string), file = file)

        for run, run_config in run_configs.items():

            times.append(timeit.timeit(
                stmt="stmt_function(**method_config, **run_config)",
                setup="setup_function(**method_config, **run_config)",
                number=iterations,
                globals={"stmt_function": stmt_function,
                         "setup_function": setup_function,
                         "method_config": method_config,
                         "run_config": run_config}
            ) / iterations)

            run_string = f"Run: {run} ({run_counter}/{num_run}) Time: {times[-1]:.6f}"
            if verbose:
                print(run_string)
            if file is not None:
                print(run_string, file=file)

            if times[-1] > cutoff_time:
                break
            run_counter += 1

        times_dict[method] = times
        method_counter += 1

        if verbose:
            print("-" * len(method_string), end="\n\n\n")
        if file is not None:
            print("-" * len(method_string), end="\n\n\n", file=file)

    return times_dict


# Plot configurations
fig, axs = plt.subplots(ncols=4, nrows=3,
                        layout="constrained", figsize=(30, 20))
fig.suptitle("pyCHomP Benchmarks")

# Report file configuration
if report_filename is not None:
    os.makedirs("examples/reports", exist_ok=True)
    path = f"examples/reports/{report_filename}"
    if os.path.exists(path):
        os.remove(path)
    f = open(path, "x")
else:
    f = None


"""
Cubical Sn Homology
-------------------
Computes the connection matrix of Sn embedded in `k`-dimensional cubical
space, with the number of boxes in each dimension given by `length`.
The parameter `length` must be at least two to properly embed Sn.

The experiment is run for all Sn where `n` < `k`. Each Sn is made up of
`n`-cells and their closure.
"""

def CubicalSnHomology_instantiate(n, k, length, **kwargs):
    X = CubicalComplex([length] * k)

    def grading(cell):
        if X.cell_dim(cell) > n:
            return 1
        barycenter = X.barycenter(cell)
        if barycenter[n+1:k] != [0] * (k - n - 1):
            return 1
        for d in range(n+1):
            if barycenter[d] > 2:
                return 1
        return 0

    include = set()
    for cell in X(n):
        if grading(cell) == 0: include.add(cell)
    include_grading = inclusion_grading(X, include)

    global gradX
    gradX = GradedComplex(X, include_grading)

def CubicalSnHomology_run(n, k, truncate=False, limit=False, **kwargs):
    cm = ConnectionMatrix(gradX, match_dim = (n + 1 if limit else -1),
                          truncate=truncate, max_grade=0, verbose=progress)
    result = [0] * (n + 2) if limit else [0] * (k + 1)
    result[0] = 1
    result[n] = 1
    assert cm.count()[0] == result

def CubicalSnHomology(n, k, length, truncate, limit, **kwargs):
    CubicalSnHomology_instantiate(n, k, length)
    CubicalSnHomology_run(n, k, truncate=truncate, limit=limit)

def CubicalSnHomology_cells(n, k, length, limit, **kwargs):
    X = CubicalComplex([length] * k)

    def grading(cell):
        if X.cell_dim(cell) > n:
            return 1
        barycenter = X.barycenter(cell)
        if barycenter[n+1:k] != [0] * (k - n - 1):
            return 1
        for d in range(n+1):
            if barycenter[d] > 2:
                return 1
        return 0

    total = 0
    Sn = 0
    for cell in X:
        if limit:
            if X.cell_dim(cell) > n + 1:
                continue
        if X.rightfringe(cell):
            continue
        if grading(cell) == 0:
            Sn += 1
        total += 1
    return (Sn, total)


k = 9
length = 3

CubicalSnHomology_method_configs = {
    "Default": {
        "truncate": False,
        "limit": False
    },
    "Truncated": {
        "truncate": True,
        "limit": False
    },
    "Dimension-Limited": {
        "truncate": False,
        "limit": True
    },
    "Dimension-Limited and Truncated": {
        "truncate": True,
        "limit": True
    }
}

CubicalSnHomology_run_configs = {
    str(n): {
        "n": n,
        "k": k,
        "length": length
    }
    for n in range(1, k)
}

CubicalSnHomology_times = benchmark(
    "Cubical Sn Homology",
    iterations,
    method_configs=CubicalSnHomology_method_configs,
    run_configs=CubicalSnHomology_run_configs,
    stmt_function=CubicalSnHomology,
    verbose=verbose,
    file=f
)

for method, times in CubicalSnHomology_times.items():
    axs[0, 0].plot(range(1, len(times)+1), times, label = method)
axs[0,0].set_title(f"Cubical Sn Homology (k={k}, l={length})")
axs[0,0].set_ylabel("Runtime (seconds)")
axs[0,0].set_xlabel("n")
axs[0,0].set_xticks(range(1, k))
axs[0,0].legend()

"""
Cubical Sn Homology (No Setup)
------------------------------
As above, but the timing only runs on the homology computation, not on
the instantiation of the cubical complex.
"""

CubicalSnHomologyNS_times = benchmark(
    "Cubical Sn Homology (no setup)",
    iterations,
    method_configs=CubicalSnHomology_method_configs,
    run_configs=CubicalSnHomology_run_configs,
    stmt_function=CubicalSnHomology_run,
    setup_function=CubicalSnHomology_instantiate,
    verbose=verbose,
    file=f
)

for method, times in CubicalSnHomologyNS_times.items():
    axs[0,1].plot(range(1, len(times)+1), times, label = method)
axs[0,1].set_title(f"Cubical Sn Homology (no setup, k={k}, l={length})")
axs[0,1].set_ylabel("Runtime (seconds)")
axs[0,1].set_xlabel("n")
axs[0,1].set_xticks(range(1, k))
axs[0,1].legend()


# Calculate cell count and sparsity
for method, method_config in CubicalSnHomology_method_configs.items():
    if method in ("Truncated", "Dimension-Limited and Truncated"):
        continue
    Sn_cells = []
    sparsity = []
    for run, run_config in CubicalSnHomology_run_configs.items():
        Sn, total = CubicalSnHomology_cells(**method_config, **run_config)
        Sn_cells.append(Sn)
        sparsity.append(100 * (1 - Sn/total))
    axs[0,2].plot(range(1, k), Sn_cells, label=method)
    axs[0,3].plot(range(1, k), sparsity, label=method)
axs[0,2].set_title("Sn Cell Count")
axs[0,3].set_title("Sn Sparsity")
axs[0,2].set_ylabel("Cells")
axs[0,3].set_ylabel("Sparsity (%)")
axs[0,2].set_xlabel("n")
axs[0,3].set_xlabel("n")
axs[0,2].set_xticks(range(1,k))
axs[0,3].set_xticks(range(1,k))
axs[0,2].legend()
axs[0,3].legend()


"""
Cubical S1 Homology
-------------------
Computes the connection matrix of S1 in cubical spaces of varying
dimensions, with the number of boxes in each dimension given by
`length`. The parameter `length` must be at least two to properly embed
S1.

The experiment is run for all dimensions `k` where `k` < `k_max`. S1 is
made up of `1`-cells and `0`-cells only.
"""

k_max = 11
length = 3

CubicalS1Homology_method_configs = {
    "Default": {
        "truncate": False,
        "limit": False
    },
    "Truncated": {
        "truncate": True,
        "limit": False
    },
    "Dimension-Limited": {
        "truncate": False,
        "limit": True
    },
    "Dimension-Limited and Truncated": {
        "truncate": True,
        "limit": True
    }
}

CubicalS1Homology_run_configs = {
    str(dim): {
        "n": 1,
        "k": dim,
        "length": length
    }
    for dim in range(2, k_max)
}

CubicalS1Homology_times = benchmark(
    "Cubical S1 Homology",
    iterations,
    method_configs=CubicalS1Homology_method_configs,
    run_configs=CubicalS1Homology_run_configs,
    stmt_function=CubicalSnHomology,
    verbose=verbose,
    file=f
)

for method, times in CubicalS1Homology_times.items():
    axs[1,0].plot(range(2, len(times)+2), times, label = method)
axs[1,0].set_title(f"Cubical S1 Homology (l={length})")
axs[1,0].set_ylabel("Runtime (seconds)")
axs[1,0].set_xlabel("k")
axs[1,0].set_xticks(range(2, k_max))
axs[1,0].legend()


"""
Cubical S1 Homology (No Setup)
------------------------------
As above, but the timing only runs on the homology computation, not on
the instantiation of the cubical complex.
"""

CubicalS1HomologyNS_times = benchmark(
    "Cubical S1 Homology (no setup)",
    iterations,
    method_configs=CubicalS1Homology_method_configs,
    run_configs=CubicalS1Homology_run_configs,
    stmt_function=CubicalSnHomology_run,
    setup_function=CubicalSnHomology_instantiate,
    verbose=verbose,
    file=f
)

for method, times in CubicalS1HomologyNS_times.items():
    axs[1,1].plot(range(2, len(times)+2), times, label = method)
axs[1,1].set_title(f"Cubical S1 Homology (no setup, l={length})")
axs[1,1].set_ylabel("Runtime (seconds)")
axs[1,1].set_xlabel("k")
axs[1,1].set_xticks(range(2, k_max))
axs[1,1].legend()


# Calculate cell count and sparsity
for method, method_config in CubicalS1Homology_method_configs.items():
    if method in ("Truncated", "Dimension-Limited and Truncated"):
        continue
    S1_cells = []
    sparsity = []
    for run, run_config in CubicalS1Homology_run_configs.items():
        Sn, total = CubicalSnHomology_cells(**method_config, **run_config)
        S1_cells.append(Sn)
        sparsity.append(100 * (1 - Sn/total))
    axs[1,2].plot(range(2, k_max), S1_cells, label=method)
    axs[1,3].plot(range(2, k_max), sparsity, label=method)
axs[1,2].set_title("S1 Cell Count")
axs[1,3].set_title("S1 Sparsity")
axs[1,2].set_ylabel("Cells")
axs[1,3].set_ylabel("Sparsity (%)")
axs[1,2].set_xlabel("k")
axs[1,3].set_xlabel("k")
axs[1,2].set_xticks(range(2,k_max))
axs[1,3].set_xticks(range(2,k_max))
axs[1,2].legend()
axs[1,3].legend()


"""
Cubical Top S1 Homology
------------------------
Computes the connection matrix of S1 in cubical spaces of varying
dimensions, with the number of boxes in each dimension given by
`length`. The parameter `length` must be at least four to properly embed
the top-dimensional S1.

The experiment is run for all dimensions `k` where `k` < `k_max`. S1 is
made up of top-dimensional cells except the one at coordinates
(1, 1, ..., 1).
"""

def FullCubicalS1Homology_instantiate(k, length, **kwargs):
    X = CubicalComplex([length] * k)

    def top_grading(cell):
        if X.coordinates(cell)[2:k] != [0] * (k - 2):
            return 1
        if X.coordinates(cell)[0:2] == [1, 1]:
            return 1
        return 0

    grading = construct_grading(X, top_grading)

    global gradX
    gradX = GradedComplex(X, grading)

def FullCubicalS1Homology_run(k, truncate=False, **kwargs):
    cm = ConnectionMatrix(gradX, truncate=truncate, max_grade=0, verbose=progress)
    result = [0] * (k + 1)
    result[0] = 1
    result[1] = 1
    assert cm.count()[0] == result

def FullCubicalS1Homology(k, length, truncate, **kwargs):
    FullCubicalS1Homology_instantiate(k, length)
    FullCubicalS1Homology_run(k, truncate=truncate)

def FullCubicalS1Homology_cells(k, length, **kwargs):
    X = CubicalComplex([length] * k)

    def top_grading(cell):
        if X.coordinates(cell) == [1] * k:
            return 1
        return 0

    grading = construct_grading(X, top_grading)

    total = 0
    S1_cells = 0
    for cell in X:
        if X.rightfringe(cell): continue
        if grading(cell) == 0: S1_cells += 1
        total += 1
    return (S1_cells, total)


k_max = 9
length = 4

FullCubicalS1Homology_method_configs = {
    "Default": {
        "truncate": False
    },
    "Truncated": {
        "truncate": True
    }
}

FullCubicalS1Homology_run_configs = {
    str(dim): {
        "k": dim,
        "length": length
    }
    for dim in range(2, k_max)
}

FullCubicalS1Homology_times = benchmark(
    "Cubical Top S1 Homology",
    iterations,
    method_configs=FullCubicalS1Homology_method_configs,
    run_configs=FullCubicalS1Homology_run_configs,
    stmt_function=FullCubicalS1Homology,
    verbose=verbose,
    file=f
)

for method, times in FullCubicalS1Homology_times.items():
    axs[2,0].plot(range(2, len(times)+2), times, label = method)
axs[2,0].set_title(f"Cubical Top S1 Homology (l={length})")
axs[2,0].set_ylabel("Runtime (seconds)")
axs[2,0].set_xlabel("k")
axs[2,0].set_xticks(range(2, k_max))
axs[2,0].legend()


"""
Cubical Top S1 Homology (No Setup)
------------------------------
As above, but the timing only runs on the homology computation, not on
the instantiation of the cubical complex.
"""

FullCubicalS1HomologyNS_times = benchmark(
    "Cubical Top S1 Homology (no setup)",
    iterations,
    method_configs=FullCubicalS1Homology_method_configs,
    run_configs=FullCubicalS1Homology_run_configs,
    stmt_function=FullCubicalS1Homology_run,
    setup_function=FullCubicalS1Homology_instantiate,
    verbose=verbose,
    file=f
)

for method, times in FullCubicalS1HomologyNS_times.items():
    axs[2,1].plot(range(2, len(times)+2), times, label = method)
axs[2,1].set_title(f"Cubical Top S1 Homology (no setup, l={length})")
axs[2,1].set_ylabel("Runtime (seconds)")
axs[2,1].set_xlabel("k")
axs[2,1].set_xticks(range(2, k_max))
axs[2,1].legend()


# Calculate cell count and sparsity
for method, method_config in FullCubicalS1Homology_method_configs.items():
    if method in ("Truncated",): continue
    S1_cells = []
    sparsity = []
    for run, run_config in FullCubicalS1Homology_run_configs.items():
        Sn, total = FullCubicalS1Homology_cells(**method_config, **run_config)
        S1_cells.append(Sn)
        sparsity.append(100 * (1 - Sn/total))
    axs[2,2].plot(range(2, k_max), S1_cells, label=method)
    axs[2,3].plot(range(2, k_max), sparsity, label=method)
axs[2,2].set_title("Top S1 Cell Count")
axs[2,3].set_title("Top S1 Sparsity")
axs[2,2].set_ylabel("Cells")
axs[2,3].set_ylabel("Sparsity (%)")
axs[2,2].set_xlabel("k")
axs[2,3].set_xlabel("k")
axs[2,2].set_xticks(range(2,k_max))
axs[2,3].set_xticks(range(2,k_max))
axs[2,2].legend()
axs[2,3].legend()


# Save plots
if plot_filename is not None:
    os.makedirs("examples/plots", exist_ok=True)
    fig.savefig(f"examples/plots/{plot_filename}")

if f is not None:
    f.close()
