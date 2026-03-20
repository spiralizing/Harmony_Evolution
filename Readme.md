### This repository contains the code for the paper: 
## "Harmony evolution and innovation in western classical music"

### Contents
```plaintext
├── LICENSE
├── Notebooks
│   ├── ChordAccuracy.ipynb: Includes the code used to compute the chord accuracy for both methods CoE and KS
│   ├── KeyAccuracy.ipynb:  Includes the code used to compute the Key accuracy for both methods CoE and KS
│   ├── KeyandNovelty.ipynb: Code to compute keys, uncertainties, diversities and novelties in the study.
│   ├── KeyFindingExample.ipynb: Example of how the CoE algorithm works
│   └── Plots.ipynb: Code used to generate the plots for the main content in the paper.
└── Readme.md
```

### Non-registry dependency (make sure to install this first):

Before using the Notebooks, make sure you have the package `MusicSpiralRepresentation.jl` in your Julia environment.

To install it you can run on your Julia REPL:
```bash
pkg> add https://github.com/spiralizing/MusicSpiralRepresentation.jl
```

