# PolyMC

### Monte Carlo code for single molecule DNA simulations

Generates equilibrium configrations within a variety of different canonical ensembles (e.g. fixed linking number, constant stretching force).

----
----
## Installation


### Linux
----

Compiling PolyMC requires make and a few c++ libraries.

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

Running a simulation requires an **input file**, an interaction **database file (IDB)** and a **sequence file**. All simulation details, such as the type of simulation (mode) and the number of Monte Carlo steps, are specified in the input file. Alternatively, all arguments specified in the input file may also be given via command line. Command line arguments will overwrite arguments given in the input file. 

The input file has to be passed as a commandline argument via the flag -in.

The minimal command to run PolyMC is
```console
./PolyMC -in input_file
```

Example input, IDB and sequence files are provided in the directory RunScripts.

----
## Input File and commandline arguments

Most arguments may be passed either via the command line or via the input file. Exceptions will be discussed below. 
In the input file, arguments have to be assigned as
```language
argument_name = value
```
Via command line, the argument name has to be given as a flag
```console
-argument_name value
```

If an arguments are provided both on via command line and input file **command line arguments take precedence**.

----
### Simulation Setup


----
### Interactions

- IDB: 

    Filename of provided interaction database file. For more information on IDB files see below

- sequence:

    Provided sequence file. See below.

- EV: (default: 0)
    
    Excluded volume diameter of the chain. Deactivated if set to zero. For more information see below.

- EV_check_crossings (default: 1)

    If true (1) moves that requires chain crossings are disallowed. This allows for fixed linking number simulations.

- subtract_T0 (default: 0)

    If set to 1, excess rotational strains are defined as the deviation of the euler vector from the static vector defined by the sequence dependent argument vec in the IDB file. Setting this argument to zero will account for static components in the rotation group prior to defining the euler vector. (for more info see [Ref 5](#es_phd))





----
### Simulation Parameters


----
### Output


----
### Dumps






- steps: (default: 0)

    Number of Monte Carlo step iterations.

- equi: (default: 0)

    Number of Monte Carlo equilibration steps. During these steps the output is disabled.

- num_bp: (default: 0)

    Number of monomers

- sigma: (default: 0)
    
    Supercoiling density 
    $$\sigma = \Delta Lk / Lk_0 = (Lk - Lk_0) / Lk_0$$ 
of the intial configuration. Lk_0 is determined by the specified intrinsic twist (via the vec argument in the IDB file). If no static twist is defined intrinsic twist is defined by helical_repeat (see next).

- helical_repeat (default: 3.57)

    Helical repeat length (in nm). Used for initiating configurations with intensive topological strain sigma. 

- force: (default: 0)

    Linear stretching force (in units of pN). Will be ignored if a closed molecule is considered

- torque: (default: 0)

    Torque (in units of pN.nm) applied to the terminal triad. 

- T: (default: 300)
    
    Temperature in Kelvin


- seed: (default: -1)

    Seed for random number generator. If set to -1 a random seed will be generated.

- dump_dir:

    Path and basename for output files. All output files will share this basename complemented with the appropriate extension. Make sure that the specified directory exists as **the directory will not be created**.

- copy_input (default: 1)

    If set to 1 all input files are included in the outfiles specified by dump_dir. The input file will be generated based on arguments provided in input file and command line. 

- 3.57

 
----
## Excluded Volumes

If EV is set to zero no excluded volume will be considered. Excluded volumes are implemented via hard sphere repulsion between beads of diameter EV. The number of excluded volume beads does not necessarily equal the 

IF EV is smaller or equal to EV/2 not every monomers will be associated with an excluded 

----
## Interaction Database (IDB)


----
## Sequence File


----
## Electrostatics


----
## Holonomic constraints
Constraints used in [4].



## Publications

1. E. Skoruppa, A. Voorspoels, J. Vreede, and E. Carlon. [Length-scale-dependent elasticity in DNA from coarse-grained and all-atom models](https://doi.org/10.1103/PhysRevE.103.042408).
*Phys. Rev. E*, 103:042408, 2021

2. M. Segers, E. Skoruppa, J. A. Stevens, M. Vangilbergen, A. Voorspoels, and E. Carlon. [Comment on “Flexibility of short DNA helices with finite-length effect: From base pairs to tens of base pairs” [J. Chem. Phys. 142, 125103 (2015)]](https://doi.org./10.1063/5.0055349). *J. Chem. Phys.*, 155(2):027101, 2021

3. E. Skoruppa and E. Carlon. [Equilibrium fluctuations of DNA plectonemes](https://doi.org/10.1103/PhysRevE.106.024412). *Phys. Rev. E*, 106:024412, 2022

4. W. Vanderlinden, E. Skoruppa , P. Kolbeck, E. Carlon, and J. Lipfert. [DNA fluctuations reveal the size and dynamics of topological domains](https://doi.org/10.1093/pnasnexus/pgac268). *PNAS Nexus*, 1:pgac268, 2022

5. <a name="es_phd"></a>E. Skoruppa, E. Carlon [Physical Modeling of DNA and DNA-Protein Interactions](https://kuleuven.limo.libis.be/discovery/fulldisplay?docid=lirias3955698&context=SearchWebhook&vid=32KUL_KUL:Lirias&lang=en&search_scope=lirias_profile&adaptor=SearchWebhook&tab=LIRIAS&query=creator%2Cexact%2CU0118787%2CAND&facet=creator%2Cexact%2CU0118787&mode=advanced). PhD thesis, KU Leuven, 2020

 



