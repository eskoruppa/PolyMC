# PolyMC

### Monte Carlo code for single molecule DNA simulations

----
----
## Installation


### Linux
----

Compiling PolymC requires make and a few c++ libraries.

Install c++ libraries 
```console
sudo apt install -y liblapack-dev
sudo apt install -y libblas-dev
sudo apt install -y libboost-dev
sudo apt install -y libarmadillo-dev
```

Compile via make
```console
make -f MakePolyMC
```

To compile using N cores run 
```console
make -f MakePolyMC -j N
```


The executable can be found in 
```console
./bin/Release/PolyMC
```


----
----
## How to run

Running a simulation requires an **input file**, a interation **database file (IDB)** and a **sequence file**. All simulation details, such as the type of simulation (mode) and the number of Monte Carlo steps are specified in the input file. Alternatively, all arguments specified in the input file may also be given via commandline. Commandline arguments will overwrite arguments given in the input file. 

The input file has to be passed as a commandline argument via the flag -in.

The minimal command to run PolyMC is
```console
./PolyMC -in input_file
```

Example input, IDB and sequence files are provided in the directory RunScripts.

----
----
## Input File and commandline arguments

Most arguments may be passed either via the commandline or via the input file. Exception will be discussed below. 
In the input file arguments have to be assigned as
```language
argument_name = value
```
Via commandline the argument name has to be given as a flag
```console
-argument_name value
```

### Most important arguments




## Publications

1. E. Skoruppa, A. Voorspoels, J. Vreede, and E. Carlon. [Length-scale-dependent elasticity in DNA from coarse-grained and all-atom models](https://doi.org/10.1103/PhysRevE.103.042408).
*Phys. Rev. E*, 103:042408, 2021

2. M. Segers, E. Skoruppa, J. A. Stevens, M. Vangilbergen, A. Voorspoels, and E. Carlon. [Comment on “Flexibility of short DNA helices with finite-length effect: From base pairs to tens of base pairs” [J. Chem. Phys. 142, 125103 (2015)]](https://doi.org./10.1063/5.0055349). *J. Chem. Phys.*, 155(2):027101, 2021

3. E. Skoruppa and E. Carlon. [Equilibrium fluctuations of DNA plectonemes](https://doi.org/10.1103/PhysRevE.106.024412). *Phys. Rev. E*, 106:024412, 2022

4. W. Vanderlinden , E. Skoruppa , P. Kolbeck, E. Carlon, and J. Lipfert. [DNA fluctuations reveal the size and dynamics of topological domains](https://doi.org/10.1093/pnasnexus/pgac268). *PNAS Nexus*, 1:pgac268, 2022

 

