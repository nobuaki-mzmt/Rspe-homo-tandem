# README

This repository provides access to the tracking data and analysis code used for the manuscript

#### Ancestral sex-role plasticity facilitates the evolution of same-sex sexual behaviour

by Nobuaki Mizumoto<sup>1</sup>, Thomas Bourguignon<sup>1, 2</sup>, and Nathan W. Bailey<sup>3</sup>.

<sup>1</sup> Okinawa Institute of Science & Technology Graduate University, Onna-son, Okinawa, Japan <br />
<sup>2</sup> Faculty of Tropical AgriSciences, Czech University of Life Sciences, Prague, Czech Republic <br />
<sup>3</sup> School of Biology, University of St Andrews, St Andrews, U.K. <br />

## Table of Contents
* [README](./README.md) - this file
* [data](./data) - folder containing raw data and created rda data during analysis
  * [phylogeny](./data/phylogeny) - phylogenetic tree for termites
  * [raw_movement](./data/raw_movement) - data for the main experiments (tracking results)
  * [rda](./data/rda) - storage for the rda data created during analysis, and rda version of other data (including grouptandem, nestmate xperiments, and comparative dataset of tandem)
  * [sim](./data/rda) - simulation results. created using "RsHomoSepSearch.cpp" 
  * [Grouptandem.xlsx] - data for group observation
  * [NestmateTandem.xlsx] - data for nest mate experiments
* [scripts](./scripts) - folder containing all script for the analysis
  * [Output.R](./scripts/Output.R) - Codes for outputting the results of movement analysis (figures, statistics, tables)
  * [Phylogeny.R](./scripts/Phylogeny.R) - Codes for phylogenetic analysis
  * [Preprocess.R](./scripts/Preprocess.R) - Codes for preprocessing data sets (run before Output.R)
  * [Sources.R](./scripts/Sources.R) - Codes for packages and functions
  * [SensitivityAnalysis.R](./scripts/SensitivityAnalysis.R) - Codes for sensitivity analysis
* [HomoTandemRspe.Rproj](./HomoTandemRspe.Rproj) - R project file
* [RsHomoSepSearch.cpp](./RsHomoSepSearch.cpp) - C++ codes for simulations

## Usage notes
To reproduce our analysis, open [HomoTandemRspe.Rproj](./HomoTandemRspe.Rproj) in [RStudio](https://www.rstudio.com/). 
[Preprocess.R](./scripts/Preprocess.R) needs to run before [Output.R](./scripts/Output.R) and [SensitivityAnalysis.R](./scripts/SensitivityAnalysis.R).
For each script, load all of the function first. Results will be produced in /plot/ or printed in the Console.

## Movement observations
We recorded movement patterns of termites during tandem run (heterosexual, female-female, and male-male) and after separation from the tandem (as a leader or as a follower). Movement trajectories were extracted from the videos using [UMATracker](https://ymnk13.github.io/UMATracker/), and the coordinates file were storaged in [raw_movement](./data/raw_movement).

## Simulations
Simulations were performed using parameters obtained from behavioral analysis (parameters will be outputted in [sim](./data/rda)). Simulations were implemented in Microsoft Visual Studio C++ 2019, with a script [RsHomoSepSearch.cpp](./RsHomoSepSearch.cpp). 

## Phylogenetic comparative analysis
We used the phylogenetic tree in [Mizumoto and Bourguignon 2021](https://doi.org//10.1098/rspb.2021.1458). The tandem information was obtained mainly from systematic search of the literature, following PRISMA statements (see Fig. S9). All the information can be found in Table S3 and S4.

See the Method and Supplementary Materials for more details.
