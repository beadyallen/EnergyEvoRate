# Overview

Repository of data and code for the paper "Diversity and Metabolic Energy in Bacteria".

Figures and tables from the manuscript may be reproduced using 'Analysis and figures.R', and should only require reasonably standard R bioinformatics packages installed such as `phyloseq`, `ape` and `vegan`. 

A more complete analysis including more diversity plots and metrics is included in the `Complete workflow - Alpha, Theta and Phylo diversity.ipynb` notebook which can be examined with the IRKernel and a Jupyter Notebook environment.

To run the pipeline with different input data, the `NMGS` R package will need installing, as well as the capability to run the actual [NMGS code](available at https://github.com/microbiome/NMGS).

`Metabolic energy estimation.ipynb` in the `metabolics` folder contains details of yield and metabolic energies used in the paper, along with code to automatcially balance stoichiometries for arbitrary catabolic and anabolic reactions.


