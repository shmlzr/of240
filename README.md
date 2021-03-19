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



