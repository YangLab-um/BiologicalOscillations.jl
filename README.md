# BiologicalOscillations

[![Build Status](https://github.com/ftavella/BiologicalOscillations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ftavella/BiologicalOscillations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ftavella.github.io/BiologicalOscillations.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ftavella.github.io/BiologicalOscillations.jl/dev)

A software package for researchers working with biological oscillations.

## Description

Welcome to BiologicalOscillations.jl. This package provides a set of tools for simulating, analyzing and visualizing oscillatory phenomena in biology.

Biological oscillations are ubiquitous among living organisms, playing critical roles in regulating key physiological processes, such as [cell cycles](https://morganlab.ucsf.edu/cell-cycle-principles-control), [circadian rhythms](https://nigms.nih.gov/education/fact-sheets/Pages/circadian-rhythms.aspx), and [neural function](https://en.wikipedia.org/wiki/FitzHugh%E2%80%93Nagumo_model). These oscillations are manifested through different signals such as cyclic changes in gene expression, electrical activity in neurons, and the periodic contraction and relaxation of muscle cells.

Despite the importance of biological oscillations, simulating and analyzing these type of signals can be challenging due to their inherent dynamic nature. By leveraging the [SciML Julia ecosystem](https://sciml.ai/), BiologicalOscillations.jl aims to simplify the process of driving insights from oscillatory data through its set of tools and algorithms.

The package includes tools for simulating protein interaction networks (PINs), gene regulatory networks (GRNs), and several important oscillatory models from the literature. Additionally, it provides functionality for detecting oscillations on differential equation models as well as methods for characterizing the generated solutions. 

The package relies heavily on [Catalyst.jl](https://github.com/SciML/Catalyst.jl), [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/), [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl), [DSP.jl](https://github.com/JuliaDSP/DSP.jl), and [Peaks.jl](https://github.com/halleysfifthinc/Peaks.jl) to name a few. Many of it’s functionalities wouldn’t be possible without these amazing packages. Be sure to check them out.

Overall, BiologicalOscillations.jl is a powerful tool for researchers studying biological oscillations, and it is designed to be flexible and extensible, so that users can adapt it to their specific research needs. We hope that this package will help advance our understanding of biological oscillations.

BiologicalOscillations.jl is developed and maintained by the [Yang Lab at the University of Michigan](http://www-personal.umich.edu/~qiongy/).