# README #

This repository contains the code and datasets required to reproduce the results of the paper Robust Stability Estimation in Correlated Space, accepted at ECML/PKDD2021.


* The file runSimu.R runs the simulations of secion 6.1. To launch these simulations, run the function runSim(index) with index being the index of the simulation (see paper). Results are saved in pdf files with names simu<INDEX>.pdf.
* The file runXP.R runs the experiments of section 6.2. The experiments presented in the paper take several tens of hours to execute (run paper() function). To faciliate the use of our software, the function demo() runs the different selection methods on the smallest dataset, alon, for M=5 selection runs only. This takes several minutes (~3) to execute. These two functions save the results in pdf files, respectively paper-phi.pdf and paper-phi_msi.pdf, and test-phi.pdf and test-phi_msi.pdf.

