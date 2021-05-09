# microclimc
**_Below, within, and above canopy microclimate modelling with R_** 

Climate strongly influences ecological patterns and processes at scales ranging from local to global. Studies of ecological responses to climate usually rely on data derived from weather stations, where temperature and humidity may differ substantially from that in the microenvironments in which organisms reside. To help remedy this, we present a model that leverages first principles physics to predict microclimate above, within, and below the canopy in any terrestrial location on earth, made freely available as an R software package. The model can be run in one of two modes. In the first, heat and vapour exchange within and below canopy are modelled as transient processes, thus accounting for fine temporal-resolution changes. In the second, steady-state conditions are assumed, enabling conditions at hourly intervals or longer to be estimated with greater computational efficiency. We validated both modes of the model with empirical below-canopy thermal measurements from several locations globally, resulting in hourly predictions with mean absolute error of 2.77 °C and 2.79 °C for the transient and steady-state modes respectively. Alongside the microclimate model, several functions are provided to assist data assimilation, as well as different parameterizations to capture a variety of habitats, allowing flexible application even when little is known about the study location. The model's modular design in a programming language familiar to ecological researchers provides easy access to the modelling of site-specific climate forcing, in an attempt to more closely unify the fields of micrometeorology and ecology.
    
Citation: Maclean, I. M. D., & Klinges, D. H. (2021). Microclimc: A mechanistic model of above, below and within-canopy microclimate. Ecological Modelling, 451, 109567. https://doi.org/10.1016/j.ecolmodel.2021.109567

## Installation

`#install.packages("devtools")`

`devtools::install_github("ilyamaclean/microclimc")`

This package builds upon software provided in the [`microctools`](https://github.com/ilyamaclean/microctools) and [`microclima`](https://github.com/ilyamaclean/microclima) packages. For further documentation on `microclima` please visit the [wiki](https://github.com/ilyamaclean/microclima/wiki).

## Issues?

If you have any problems with the package or wish to suggest improvements, then please open an issue [here](https://github.com/ilyamaclean/microclimc/issues).
