# of240

Installation


For the frozen solver
```
cd solver/frozenSimpleFoam
wmake 
```

For the frozen turbulence model
```
cd turbulenceModels/frozenTurbModels
wmake 
```

add this line to the controlDict
```
libs ( "libFrozenTurbModels.so" );
```

In order to run it, do the following


First change
```
constant/RASProperties
RASModel        frozenOmegaSST_modPk;

```

Than run
```
frozenSimpleFoam

```

## Relevant Publications

Schmelzer, M., Dwight, R.P. & Cinnella, P.: Discovery of Algebraic Reynolds-Stress Models Using Sparse Symbolic Regression. 
Flow Turbulence Combust 104, 579–603, 2020
https://doi.org/10.1007/s10494-019-00089-x

Ismaïl Ben Hassan Saïdi, Martin Schmelzer, Paola Cinnella, Francesco Grasso: CFD-driven symbolic identification of algebraic Reynolds-stress models,
Journal of Computational Physics, Volume 457, 2022
https://doi.org/10.1016/j.jcp.2022.111037

Further applications:
https://www.researchgate.net/profile/Martin-Schmelzer


##
2021, Martin Schmelzer, m.schmelzer@tudelft.nl

