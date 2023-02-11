# Dynamic Smagorinsky (DynSmag)
The Dynamin Smagorinsky SGS model with local averaging, based onThe Dynamin Smagorinsky SGS model with local averaging, based on the following works:

    1- [M. Germano, U. Piomelli, P. Moin, and W. H. Cabot, “A dynamic subgrid-scale eddy viscosity model,” Phys. Fluids A 3, 1760–1765 (1991).](https://aip-scitation-org.qe2a-proxy.mun.ca/doi/10.1063/1.857955)
    2- D. K. Lilly, “A proposed modification of the Germano subgrid-scale closure method,” Phys. Fluids A 4, 633–635 (1992).
 
This code had been developed for OpenFOAM version 4.1 by Ehsan Asgari and an older one by Alberto Passalacqua. This version is developed for OpenFOAM-10.

## How to use
In your OpenFOAM case add the following code line to the "system/controlDict" file:

```
libs ("/path/to/DynSmag.so");
```

and you are good to run the case. Happy computing!
