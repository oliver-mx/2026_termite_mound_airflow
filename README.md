# Learning from modelling and simulation of air flow dynamics inside termite mounds in view of low-energy houses

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/1359170631.svg)](https://doi.org/10.5281/zenodo.22708206)


This repository contains information and code to reproduce the results presented in the article
```bibtex
@online{marx26,
  title={Modelling and simulation of air flow dynamics inside termite mounds},
  author={Marx, Oliver P and Gasser, Ingenuin and Schmidgall, Annika},
  year={2026},
  howpublished={\url{https://github.com/oliver-mx/2026_termite_mound_airflow}},
  version={0.1.1},
  doi = {https://doi.org/10.5281/zenodo.22708206}
}
```


## Abstract

Termite mounds arecomplex animal-built structures with varied architectural designs depending on the species and respective habitats.
One key function is to provide ventilation for the underground nest. 
This regulates the temperature of the nest and diffuses the metabolic gases produced inside. 
A one-dimensional mathematical model is presented that describes the airflow inside of a termite mound. 
This thermo fluid dynamic model originates from the Euler equations of gas dynamics in a low Mach number regime. 
Numerical simulations are performed to study full 24 hour day-night air flow cycles and the occurrence of day-night air flow cycles for differently sized termite mounds. 
An optimal mound geometry is determined, which achieves stable nest temperature control. 
The results are in good agreement with measured data from the literature.


## Numerical experiments

In order to reproduce the numerical experiments presented in this article, you need to install [Julia](https://julialang.org/). 
The numerical experiments presented in this article were performed using Julia v1.12.6.
Download this repository, e.g., by cloning it with `git` or by downloading an archive via the GitHub interface.
Then, start Julia in the `code` directory of this repository and follow the instructions described in the `README.md` file therein.


## Authors

- Oliver P Marx (University of Hamburg, Germany)
- Ingenuin Gasser (University of Hamburg, Germany)
- Annika Schmidgall (University of Hamburg, Germany)


## License

The code in this repository is published under the MIT license, see the `LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!