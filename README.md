# Simulation of marRAB *e. coli* gene circuit

[![DOI](https://zenodo.org/badge/9321/spencermamer/marRAB-SSA.svg)](http://dx.doi.org/10.5281/zenodo.14307)

This is an implementation of the Gillespie stochastic simulation algorithm (**SSA**) as applied to the marRAB-rob gene circuit, which I put together a year or so ago. The code contains the kernel of a general SSA toolset in Java (class-based definining of reactions, their species and rate parameters, and so on. It is _not_ intended for those purposes, because I was not yet attempting to develop a stochastic reaction modeling toolset. I was attempting, as a learning exercise, to reproduce the pulsing behavior modeled in a paper published by Garcia-Bernardo et al in 2013, entitled "Tunable Stochastic Pulsing in the Escherichia coli Multiple Antibiotic Resistance Network from Interlinked Positive and Negative Feedback Loops.' (DOI: 0.1371/journal.pcbi.1003229)
