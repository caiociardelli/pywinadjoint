# PyWinAdjoint

`PyWinAdjoint` is a set of scripts for carrying out time-window selection and computing exponentiated phase adjoint sources for full-waveform inversion.

Author: Caio Ciardelli

If you use PyWinAdjoint, please, cite the following paper:

Ciardelli, C., Assumpção, M, Bozdağ, E., and  van der Lee, S., 2021. Adjoint Tomography of South America based on 3D Spectral-Element Seismic Wave Simulations. JGR: Solid Earth, submitted.

[![DOI](https://zenodo.org/badge/394311067.svg)](https://zenodo.org/badge/latestdoi/394311067)

## Installation

`PyWinAdjoint` requires no installation, but some requirements must be met:

* Python3.8 or later
* Obspy1.1.0 or later
* FFTW3 or later

To compile the code, just run the Makefile:

```shell
make
```

## Development

Development is hosted on GitHub in the [caio.ciardelli/pywinadjoint repository](https://github.com/caiociardelli/pywinadjoint).

## Other References
-----------------------

Please, also cite:

Yuan, Y.O., Bozdağ, E., Ciardelli, C., Gao, F., Simons, F.J., 2020. The exponentiated phase measurement, and objective-function hybridization for adjoint waveform tomography. Geophysical Journal International, v. 221, no. 2, p. 1145–1164. doi: 10.1093/gji/ggaa063,

for the exponentiated phase measurement and:

Ruan, Y., Lei, W., Modrak, R., Örsvuran, R., Bozdağ E., Tromp, J., 2019. Balancing unevenly distributed data in seismic tomography: a global adjoint tomography example. Geophysical Journal International, v. 219, no. 2, p. 1225–1236, doi:10.1093/gji/ggz356,

for the receiver weights to balance uneven station distribution.

## Usage

This package is still under development, hence it still lacks a detailed documentation. You can set the parameters in the "Config.cfg" file.

To execute the example, just run:

```shell
./run_example.bash
```

In case you set "save_binary" to "True" in the "Config.cfg" file, you need to convert the adjoint sources from binary to the ASCII format. To do that, run:

```shell
./utils/expand.bash 1 4

## Contact
-----------------------

If you have any questions, suggestions, and bug reports, you can email *caio.ciardelli@gmail.com*

