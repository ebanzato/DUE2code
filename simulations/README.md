# Simulation study

This folder contains the R code for the simulation study, which evaluates the
two-stage procedure implemented in the R package
[DUE2](https://github.com/ebanzato/DUE2).


## Folder content

| File | Description |
|---|---|
| `simulations.R` | Main script: builds the graph and the parameters, runs the simulations for k = 2, 3, 4, computes the performance measures and produces the figures. |
| `sim_functions.R` | Functions to simulate data and to compute the performance measures (`OFDR`, `power.1stage`, `fwer.2stage`, `power.2stage`). |




## How to run

Set the paths in `simulations.R`:

   - `source('path/to/sim_functions.R')`: the simulation functions;
   
   - `path/to/results/`: the folder where the results are saved and loaded
     (`TS_tpgm_2.RData`, `TS_tpgm_3.RData`, `TS_tpgm_4.RData`).

**Running time.** On a
standard PC with 7 cores, each value of k takes about 2–3 days, so the whole
simulation study takes about one week.



