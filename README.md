# BPASS_GridFit
### Code to fit measured line fluxes to a grid of photoionization models

The code aims to calculate the best-fitting ISM parameters for a given target based on its optical emission lines when compΩred to a grid of theoretical flux predictions.

The theoretical emission line grid should sample different values of metallicity $Z$, ionization parameter $\log(U)$, metal-to-dust ratio $\xi$, hydrogen density $\log(n_\text{H})$, and C/O fraction.

The code is optimized for the following emission line grid:

- Stellar population models

  BPASS

  Eldridge+2017

  https://arxiv.org/abs/1710.02154

- Photoionization code

  cloudy

  Chatzikos+2023

  https://arxiv.org/pdf/2308.06396
