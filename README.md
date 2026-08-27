# How demography shapes linkage disequilibrium with or without recombination

This repository contains the C code used for the simulation study in:

> Kerdoncuff, E., Ianni-Ravn, M. K., Roze, D., Novembre, J., Lambert, A. & Achaz, G.
> *How demography shapes linkage disequilibrium with or without recombination.*
> bioRxiv, 2025.
> https://doi.org/10.1101/2023.12.13.571342

---

## Overview

Linkage disequilibrium describes the statistical association between alleles at different loci. Its observed patterns can be affected by several evolutionary processes, including recombination, genetic drift, and demographic history.

This repository implements the simulations used to investigate these effects.

In particular, the simulations consider:

* **No recombination** between the two loci
* **A single recombination event**
* **Independent loci**, corresponding to an effectively infinite number of recombination events
* **Changes in population size** under a two-epoch demographic model
* Different measures of LD:

  * \(D\)
  * \(|D'|\)
  * \(r^2\)
  * 2-site Configuration Probabilities (CPs)

The simulations are based on genealogical trees rather than explicitly generating polymorphic sequences. 

---

## Repository structure

```text
LD_study/
├── D_etude_lib.c
├── D_etude_lib.h
├── D_stats_study_double_treeweight.c
├── run_D_stats_study_CP_treeweight.sh
└── README.md
```

### `D_etude_lib.c`

Core simulation library implementing the genealogical/coalescent machinery and functions used by the simulation programs.

### `D_etude_lib.h`

Header file containing the declarations required by `D_etude_lib.c` and the simulation programs.

### `D_stats_study_double_treeweight.c`

Main simulation program used to generate LD statistics under different demographic and recombination parameters.

### `run_D_stats_study_CP_treeweight.sh`

Shell script for compiling the simulation and running it over a range of demographic parameters.

---

## Requirements

The code is written in C and is primarily intended for macOS.

It uses platform-specific memory-allocation functions:

- On macOS, `malloc_size()` from `<malloc/malloc.h>` is used.
- On Linux, the corresponding function is `malloc_usable_size()` from `<malloc.h>`.

The current implementation therefore requires a small modification to compile on Linux.

The simulations also require:

- A C compiler (`gcc` or `clang`)
- A POSIX-compatible shell
- `bc`

### macOS

The code can be compiled directly with:

```bash
gcc D_stats_study_double_treeweight.c D_etude_lib.c \
    -o D_etude_study_treeweight.exe
```

---

## Running the simulations

The main executable accepts a set of command-line parameters controlling the sample size, demographic history, output, simulation mode, random seed, and number of iterations.

The general syntax is:

```bash
./D_etude_study_treeweight.exe \
    -N <N> \
    -T <TCHANG> \
    -C <CONTRAC> \
    -A <TREE> \
    -S <SEED> \
    -D <STATS_D> \
    -n <ITERATIONS> \
    -R <RECOMBINATION_REPORT>
```

The available parameters are:

| Option | Parameter              | Description                                                                                  |
| ------ | ---------------------- | -------------------------------------------------------------------------------------------- |
| `-N`   | `N`                    | Number of leaves (sample size) in the simulated genealogical tree.                           |
| `-T`   | `TCHANG`               | Time of the demographic change (contraction/expansion) in coalescent time units.             |
| `-C`   | `CONTRAC`              | Strength of the demographic contraction/expansion.                                           |
| `-A`   | `TREE`                 | Controls whether the simulated Newick tree is printed: `0` = do not print, `1` = print.      |
| `-S`   | `SEED`                 | Seed for the random-number generator.                                                        |
| `-D`   | `STATS_D`              | Selects the LD/statistical analysis to perform.                                              |
| `-n`   | `n`                    | Number of simulation iterations.                                                             |
| `-R`   | `RECOMBINATION_REPORT` | Controls whether results are reported separately by recombination type: `0` = no, `1` = yes. |

### `-N`: Number of leaves

```bash
-N <N>
```

Sets the number of leaves in the simulated genealogical tree, corresponding to the number of sampled chromosomes/individuals represented in the tree.

---

### `-T`: Time of demographic change

```bash
-T <TCHANG>
```

Sets the time at which the population-size change occurs.

---

### `-C`: Strength of demographic change

```bash
-C <CONTRAC>
```

Sets the strength of the demographic contraction or expansion.

---

### `-A`: Print the Newick tree

```bash
-A <TREE>
```

Controls whether the simulated Newick tree is printed.

* `0`: do not print the tree
* `1`: print the tree

Mostly use for debugging purposes.

---

### `-S`: Random-number seed

```bash
-S <SEED>
```

Sets the seed used by the random-number generator.

For reproducible simulations, it is therefore recommended to explicitly specify a non-zero seed.

---

### `-D`: Select the statistical analysis

```bash
-D <STATS_D>
```

Selects which statistical regime is simulated.

The currently defined values include:

| Value | Meaning                                                                      |
| ----: | ---------------------------------------------------------------------------- |
|   `0` | \(D_0\): no recombination                                                    |
|   `1` | \(D_1\): one recombination                                                   |
|   `8` | \(D_\infty\): infinitely many recombination events / independent genealogies |


The selected mode determines how the two-locus genealogies are generated and which quantities are subsequently calculated.

The simulations report both classical LD statistics and configuration probabilities, including:

$$
D,\quad D',\quad r,\quad r^2
$$

and

$$
CP_0,\quad CP_1,\quad CP_2,\quad CP_3.
$$

---

### `-n`: Number of iterations

```bash
-n <ITERATIONS>
```

Sets the number of Monte Carlo simulation iterations.

Increasing the number of iterations reduces Monte Carlo sampling noise but increases the computational time.

---

### `-R`: Report by recombination type

```bash
-R <RECOMBINATION_REPORT>
```

Controls whether the results are reported separately according to recombination type.

* `0`: do not report statistics separately by recombination type
* `1`: compute and report statistics separately for each recombination type

This option is useful when the goal is to examine how the different recombination configurations contribute to the resulting LD statistics.

---

## Example

A complete simulation command could be:

```bash
./D_etude_study_treeweight.exe \
    -N 10 \
    -T 0.5 \
    -C 1 \
    -A 0 \
    -S 12345 \
    -D 8 \
    -n 10000 \
    -R 0 \
    output.txt
```

This runs 10,000 simulations with:

* 10 leaves;
* a demographic change at \(T_{\mathrm{CHANG}}=0.5\);
* a contraction parameter of `1`;
* no Newick tree output;
* random seed `12345`;
* the \(D_\infty\) statistical regime (`-D 8`);
* no separate reporting by recombination type;
* results written to `output.txt`.

The output filename is supplied as the final positional argument, after the command-line options.

---


