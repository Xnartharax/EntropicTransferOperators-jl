# Entropic Transfer Operators 

This repository contains the code to reproduce most of the results of my Bachelorthesis. We provide a generic implementation of the Entropic Transfer operator technique [1] and utitilites for visualization and computation of the spectrum. Here is an overview over the implementation files:

- Single.jl: Implementation for computing a single transfer operator and its spectrum
- Persistent.jl: Functionality to compute multiple transfer operators efficiently for varying regularization terms
- SpectrumPlots.jl: Tools to visualize complex spectra
- Tracking.jl: Tracking real eigenvalues over a scale of regularizations

We apply this implementation to 5 examples:

- The circle shift map (Circle.jl)
- The Four Legs map (FourLegs.jl)
- The Lorenz system (Lorenz.jl)
- Alanine-Dipeptide (Molecule.jl)
- The El Nino oscillation (Oceanic.jl)


[1] O. Junge, D. Matthes, B. Schmitzer, Entropic Transfer Operators, https://arxiv.org/abs/2204.04901