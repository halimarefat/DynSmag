# Dynamic Smagorinsky (DynSmag)
The Dynamin Smagorinsky SGS model with local averaging, based on Germano and Lilly publications. This code had been developed for OpenFOAM version 4.1 by Ehsan Asgari and an older one by Alberto Passalacqua. The current version has been developed for OpenFOAM-10.

## How to use
In your OpenFOAM case add the following code line to the "system/controlDict" file:

```
libs ("/path/to/DynSmag.so");
```

and you are good to run the case. Happy computing!
