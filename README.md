# README for [Moser and Yared (2021)](https://www.nber.org/papers/w27062)


## Description

This repository contains the replication code and materials for [Moser and Yared (2021)](https://www.nber.org/papers/w27062).

### Replication code

The following code files are included in the repository, in the order in which they need to be run:

- **[master.m](master.m)**: Simulate an economy subject to optimal lockdown policy with and without commitment for [Moser and Yared (2021)](https://www.nber.org/papers/w27062). MATLAB master file that can also be called from Stata in master.do -- run this first!
   - **[simulations.m](simulations.m)**: Simulate the evolution of key variables and return summary statistics. Called by [master.m](master.m).
      - **[backward_induct_comm.m](backward_induct_comm.m)**: Return the continuation value and the optimal lockdown policy from time tt to TT for a government with commitment. Called by [simulations.m](simulations.m).
      - **[backward_induct_nocomm.m](backward_induct_nocomm.m)**: Return the continuation value and the optimal lockdown policy from time tt to TT for a government without commitment. Called by [simulations.m](simulations.m).
      - **[gdp_vax.m](gdp_vax.m)**: Compute the continuation value of GDP after arrival of the vaccine. Called by [simulations.m](simulations.m).
      - **[utility_vax.m](utility_vax.m)**: Compute the continuation value of welfare after arrival of the vaccine. Called by [simulations.m](simulations.m), [backward_induct_comm.m](backward_induct_comm.m), and [backward_induct_nocomm.m](backward_induct_nocomm.m).
      - **[sir_model.m](sir_model.m)**: Return next period's health state, assuming no more new infections can occur after the arrival of the vaccine. Called by [simulations.m](simulations.m), [backward_induct_comm.m](backward_induct_comm.m), [backward_induct_nocomm.m](backward_induct_nocomm.m), [gdp_vax.m](gdp_vax.m), and [utility_vax.m](utility_vax.m).
- **[master.do](master.do)**: Create time-series and comparative-statics figures for an economy subject to optimal lockdown policy with and without commitment for [Moser and Yared (2021)](https://www.nber.org/papers/w27062). Stata master file that can also call [master.m](master.m) via MATLAB shell -- run this second!
   - **[time_series.do](time_series.do)**: Create time-series figures. Called by [master.do](master.do).
   - **[comp_stat.do](comp_stat.do)**: Create comparative-statics figures. Called by [master.do](master.do).

### Materials

The repository also includes the following data used in the analysis:

- **[IO_405_industries_RED.xlsx](IO_405_industries_RED.xlsx)**: Classification of materials used as intermediate inputs based on the Input-Output Accounts Data by the [U.S. Bureau of Economic Analysis](https://www.bea.gov/industry/input-output-accounts-data).


## System requirements

All codes were run on an iMac Pro (2017) with 3 GHz 10-Core Intel Xeon W processor and 128 GB 2666 MHz DDR4 memory running macOS Catalina Version 10.15.7. The codes were executed using MATLAB R2020b Update 3 (9.9.0.1538559) for all files in .M format and Stata/MP 15.1 for Mac (64-bit Intel) for all files in .DO format. The codes are compatible with other hardware, operating systems, and software versions.


## Computation time

The expected computation time is approximately 20 minutes on a machine with the above-referenced specifications.


## Random numbers

No random numbers are used in the analysis.


## References

[Moser, Christian and Pierre Yared. 2021. "Pandemic Lockdown: The Role of Government Commitment." Accepted at the *Review of Economic Dynamics*.](MY2021.pdf)
