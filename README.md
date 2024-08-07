# PolyMC

### Monte Carlo code for single molecule DNA simulations

Generates equilibrium configurations within a variety of different canonical ensembles (e.g. fixed linking number, constant stretching force). 

For more info see [Ref 5](#es_phd), Supplement of [Ref 4](#vand22) and Appendix of [Ref 3](#skor22). Note that in these references, most algorithmic details are omitted. More detailed information will be provided in the future. 


## Table of Contents

1. [Installation](#install)

    1.1 [Linux](#install_linux)

2. [How to run](#howtorun)

3. [Input File and command line arguments](#inputfile)

    3.1 [Simulation Setup](#simsetup)

    3.2 [Interaction Setup](#interactionsetup)

    3.3 [Simulation Parameters](#simparams) 

    3.4 [Output](#output) 

    3.5 [Dump Setup](#dumpsetup) 

4. [Dumps](#dumps)

    <!-- 4.1 [XYZ](#simsetup)
    
    4.2 [Simulation State](#dump_state)

    4.3 [Thetas](#dump_thetas)

    4.4 [End-to-End Distance](#dump_endend)

    4.5 [z-extension](#dump_zext)

    4.6 [Force Extension](#dump_fext)

    4.7 [Linking number](#dump_lk)

    4.8 [Writhe Map](#dump_writhemap)

    4.9 [Persistence Length](#dump_persistencelength)

    4.10 [Tangent-tangent Correlation Function](#dump_tantancorr)

    4.11 [Restart Snapshots](#dump_restart) -->

5. [Interaction Database (IDB)](#idb)

6. [Sequence File](#seq)

7. [IOPolyMC](#iopolymc)

8. [Citing PolyMC](#citing)

9. [Publications](#publications)

10. [References](#refs)





----
----
# Installation <a name=install></a>

### Linux <a name=install_linux></a>

Compiling PolyMC requires make and a few c++ libraries.

Install c++ libraries 
```
sudo apt install -y liblapack-dev
sudo apt install -y libblas-dev
sudo apt install -y libboost-dev
sudo apt install -y libarmadillo-dev
```

Compile via make
```
make -f MakePolyMC
```

To multithread the compilation run
```
make -f MakePolyMC -j N
```


After compiling, the executable can be found at 
```
./bin/Release/PolyMC
```


----
----
# How to run <a name=howtorun></a>

Running a simulation requires an `input file`, an `interaction database file (IDB)`, and a `sequence file`. All simulation details, such as the type of simulation (mode) and the number of Monte Carlo steps, are specified in the input file. Alternatively, all arguments specified in the input file may also be given via command line. Command line arguments will overwrite arguments given in the input file. 

The input file has to be passed as a command line argument via the flag -in.

The minimal command to run PolyMC is
```
./PolyMC -in input_file
```

Example input, IDB, and sequence files are provided in the directory RunScripts.

----
----
# Input File and command line arguments <a name=inputfile></a>

Most arguments may be passed either via the command line or via the input file. Exceptions will be discussed below. 
In the input file, arguments have to be assigned as
```language
argument_name = value
```
Via command line, the argument name has to be given as a flag
```
-argument_name value
```

If arguments are provided both via the command line and input file, **command line arguments take precedence**.

----
## Simulation Setup <a name=simsetup></a>


### Simulation Modes <a name=modes></a>

- `mode`:

    Specify simulation mode. Currently, implemented modes are:

    - `open`: Linear molecule without any linking constraints. 
    - `tweezer`: linear molecule emulating the setup of magnetic tweezers as used in Refs [3](#skor22) and [4](#vand22). Terminal tangents remain aligned, and conservation of linking number may be imposed via impenetrable planes aligned with the terminal monomers. Has four different setups specified by the argument `tweezer_setup`:
        - `free`: allows for free fluctuations of the total linking number
        - `fix_link`: fixed linking number ensemble
        - `torsional_trap`: traps the linking number in a quadratic potential to allow for the torque to be measured. Corresponds to quasi-fixed linking number ensemble.
        - `torque`: constant torque ensemble 
    - `plasmid`: circular molecule. Linking number remains conserved if excluded volumes are active and EV_check_crossings=1 (see below)
    - `umbrellaplasmid`: circular molecule with variable linking number to allow for umbrella sampling the plectoneme free energy
    - `plec`: pure plectoneme simulation used in Ref [3](#skor22). 
    - `lin2d`: linear molecule in 2 dimensions
    - `closed2d`: circular molecule in 2 dimensions


    Each mode has additional mode-specific arguments. See example files in RunScripts directory.

### Restart files <a name=restart_files></a>
Simulations may be initiated from previously generated snapshots if during these simulations restart files were dumped ([see in the section on dumps](#dump_restart)). Initiation from restart file is controlled by two arguments:

- `restart`: 

    Specifies the restart file.

- `restart_snapshot`: (default: -1 -> last snapshot in file)

    Index of chosen snapshot within restart file.


----
## Interaction Setup <a name=interactionsetup></a>

- `IDB`: 

    The filename of the provided interaction database file. For more information on IDB files, see below

- `sequence`:

    Provided sequence file. See below.

- `EV`: (default: 0)
  
    Excluded volume diameter of the chain. Deactivated if set to zero. Excluded volumes are implemented via hard sphere repulsion between beads of diameter EV (See Supporting of [Ref 4](#vand2022)). Note that the number of excluded volume beads does not necessarily equal the number of chain monomers if the discretization length is smaller than the bead diameter. To improve performance, as few excluded volume beads as possible are considered.

- `EV_check_crossings` (default: 1)

    If true (1), moves that lead to effective chain crossings are disallowed. This allows for fixed linking number simulations.

- `subtract_T0` (default: 0)

    If set to 1, excess rotational strains are defined as the deviation of the Euler vector from the static vector defined by the sequence-dependent argument 'vec' in the IDB file. Setting this argument to zero will account for static components in the rotation group prior to defining the Euler vector. (for more info see [Ref 5](#es_phd))


----
## Simulation Parameters <a name=simparams></a>

- `steps`: (default: 0)

    Number of Monte Carlo step iterations.

- `equi`: (default: 0)

    Number of Monte Carlo equilibration steps. During these steps, the output is disabled.

- `num_bp`: (default: 0)

    Number of monomers

- `sigma`: (default: 0)
  
    Supercoiling density 
    $$\sigma = \Delta Lk / Lk_0 = (Lk - Lk_0) / Lk_0$$ 
    of the initial configuration. Lk_0 is determined by the specified intrinsic twist (via the 'vec' argument in the IDB file). If no static twist is defined, the intrinsic twist is defined by helical_repeat (see next).

- `helical_repeat` (default: 3.57)

    Helical repeat length (in nm). Used for initiating configurations with intensive topological strain sigma. 

- `force`: (default: 0)

    Linear stretching force (in units of pN). Will be ignored if a closed molecule is considered

- `torque`: (default: 0)

    Torque (in units of pN.nm) applied to the terminal triad. 

- `T`: (default: 300)
  
    Temperature in Kelvin


- `seed`: (default: -1)

    Seed for random number generator. If set to -1 a random seed will be generated.


----
## Output <a name=output></a>
Configure command line output.

- `print_every` (default: 100000)

    Prints the status of the simulation every so many iteration steps.

- `print_link_info` (default: 0)

    Specifies whether writhe, twist, and linking number shall be displayed. This may slow down the simulation if printing frequently to command line.

----
## Dump Setup <a name=dumpsetup></a>

- `dump_dir`:

    Path and base name for output files. All output files will share this base name, complemented with the appropriate extension. Make sure that the specified directory exists, as **the directory will not be created**.

- `copy_input` (default: 1)

    If set to 1 all input files are included in the outfiles specified by dump_dir. The input file will be generated based on arguments provided in input file and command line. 

- `app` (default: 0)

    Specifies whether output is appended to existing files or whether existing files are overwritten.

---
# Dumps <a name=dumps></a>
PolyMC offers several different types of dumps that enable configurations and observables to be saved to a file.

To activate most dumps, set the corresponding "dump_every" parameter to a positive value. This will cause the corresponding output to be saved to a file after a certain number of Monte Carlo (MC) steps. If no additional filename is specified, the dump filename will be the one specified in the "dump_dir" parameter, with the appropriate extension.

----
## XYZ <a name=dump_xyz></a>
Dumps configuration in xyz format.

- Translation Options:
    - first:    First Triad always has coordinates 0,0,0 (default)
    - COM:      Center of Mass always has coordinates 0,0,0
    - last:     Last Triad always has coordinates 0,0,0

- Representation Options:
    - simple:   One bead per triad (default)
    - fl:       First and last monomer emphasized
    - helix:    helix representation
    - dna:      double helix representation
    - triadf:   Triad representation   	

```
dump_every:         XYZn
filename:           XYZfn
translate_option:   XYZ_translate
representation:     XYZ_repr

extension:          .xyz
```

Example for input file arguments: 
```
dump_dir        = dump/test
XYZn            = 10000
XYZ_translate   = COM
```
This command will initiate the configuration to be printed to file in xyz format every 10000 steps and the center of mass of the configuration to be shifted to the origin of the coordinate system. The resulting output files will be 
```
dump/test.xyz
```

The same can be achieved via command line
```
-dump_dir dump/test -XYZn 10000 -XYZ_translate COM -XYZfn dump/xyzout
```
This command will generate the same output as the input file command. Except in this case, we don't want to use the colloquial filename given by the dump_dir argument, but instead use the unique file name specified by -XYZfn. Generates the output file (this will not automatically add the .xyz extension)
```
dump/xyzout
```

Note that only one xyz dump can be active at a time. 

----
## Simulation State <a name=dump_state></a>
Prints the state of the simulation to file. Monomer positions are printed by default. Triads and angles may be printed optionally if the corresponding flags are set. The arguments dump_triads and dump_omegas are booling flags (valid values 0 and 1)
					
    dump_every:     Stn
    filename:       Stfn
    dump_triads:    Sttriads,Sttds
    dump_omegas:    StOmegas,StOm
    
    extension:      .state

Example (input file)
```
Stn         = 1000
dump_triads = 1
dump_omegas = 1
```
Prints state including monomer positions, triads, and angles to file every 1000 steps. Generates a file with name specified by the dump_dir flag with the extension complemented by the '.state' extension.

----
## Thetas <a name=dump_thetas></a>
Prints Euler vectors (Thetas) connecting consecutive triads to file
					
    dump_every:     Thetasn
    filename:       Thetasfn
    
    extension:      .thetas

----
## End-to-End Distance <a name=dump_endend></a>
Distance between first and last monomer
	
    output: distance
    
    dump_every:     E2En 
    filename:       E2Efn
    
    extension:      .e2e

----
## z-extension <a name=dump_zext></a>
Extension along the force direction. This dump can be used at every step without significant loss of efficiency. 
	
	output: zext
	
	dump_every: Extn,  extn,  EXTn 
	filename:   Extfn, extfn, EXTfn
	
	extension: .zext

----
## Force Extension <a name=dump_fext></a>
Calculates the force-extension statistics in the direction of the force and prints them to a file. Once the simulation is complete, a single line of output is generated. By specifying a fraction of the chain, only the middle portion of the monomers will be included in the calculation. For example, setting the "fefrac" parameter to 0.5 means only the middle half of the monomers will be used. This can help to avoid finite-size effects that can arise from boundary terms.
	
	output: force number_of_measurements z z_squared contour_length
	
	dump_every: fen,  FEn   
	filename:   fefn, FEfn
	fraction:   fefrac, FEfrac (default 1, should be 0 < frac <= 1)
	
	extension: .fe

----
## Energy <a name=dump_energy></a>
Total elastic energy in the system in units of kT
	
```
output: energy

dump_every: En 
filename:   Efn

extension: .en
```

----
## Linking number <a name=dump_lk></a>
Prints writhe and twist.

Options for writhe: 
- exact:  Langowski method 1a (see [Klenin 2000](#klen00))
- fuller: Fuller single sum (see [Fuller 1978](#full78))
- both:   Both Langowski and Fuller (3 arguments: tw lang_wr full_wr)
fsdfsdf

```
dump_every:     LKn 
filename:       LKfn
writhe option:	LKoptn
chain fraction: LKfrac

extension:      .lk
```

----
## Writhe Map <a name=dump_writhemap></a>
Prints the writhe map, the pairwise components of the double sum (see Ref [3](#skor22)) to file.

    dump_every:     WMn
    filename:       WMfn
    segment_size:   WMseg
    
    extension:      .wm

----
## Persistence Length <a name=dump_persistencelength></a>
Calculates the persistence length via the tangent-tangent correlation function using the inversion (see Ref [5](#es_phd))

$$
l_b(m) = \frac{-am}{\log\langle \hat{\mathbf{t}}_i \cdot \hat{\mathbf{t}}_ {i+m} \rangle}
$$

In this method, 'm' denotes the number of monomers by which the tangents are displaced, and 'a' represents the discretization length. The 'maxdist' argument sets the maximum value of 'm' that is considered.
					
    dump_every:     LBn 
    filename:       LBfn
    maxdist:        LBdist
    
    extension:      .lb

----
## Tangent-tangent Correlation Function <a name=dump_tantancorr></a>
Analogous to persistence length
					
    dump_every:     TCn
    filename:       TCfn
    maxdist:        TCdist
    
    extension:      .tancor

----
## Restart Snapshots <a name="dump_restart"></a>
Restart files allow the simulation to be resumed from a previously generated snapshots. 
					
    dump_every:     Restartn, restartn
    filename:       Restartfn, restartfn
    
    extension:      .restart


----
----
# Interaction Database (IDB) <a name="idb"></a>
The interaction database file specifies the Hamiltonian of the polymer chain as well as the ground state configuration or static rotational structure. Every Hamiltonian is expressed in terms of the rotational strain fields, i.e. the rotational deformations away from the ground state.

The first few arguments specify the general setup:

```
interaction_range  = 0
monomer_types      = a
discretization     = 0.34
avg_inconsist      = 0
```
The argument 'interaction_range' sets the range of the couplings of the rotational deformations. If set to zero, interactions are strictly local. In this case, the energy is fully determined by all possible dimers. For non-zero interaction range, larger oligomers have to be considered. In general, the size of the oligomers to be specified is 

$$
2\times(\mathrm{interaction} \; \mathrm{range}+1).
$$


In the remainder of the file, all possible oligomers have to be specified. The form of the energy of the diagonal and off-diagonal terms may be specified via the leading indicator. For example, 'stiffmat' indicates a 3x3 stiffness matrix. 

Example: If the interaction range is 1 and the chain consists only of type 'a' monomers, only one interaction has to be specified:

```
aaaa
    stiffmat    5  0 0 0 5  0 0 0 10
    stiffmat    40 0 0 0 40 0 0 0 100
    stiffmat    5  0 0 0 5  0 0 0 10
    vec	    	0 0 36
```

To operate this system, three stiffness matrices must be specified: The first matrix represents the coupling of the middle step to the left neighbor, the second matrix specifies the diagonal component of the middle step, and the third matrix represents the couplings of the middle step to the right neighbor.
The nine numbers following the 'stiffmat' indicator correspond to the nine entries of a 3x3 matrix, specified row-wise (i.e., filling the first, second, and third row consecutively).

If not carefully chosen, the left couplings may not match the right couplings. Such inconsistencies may be averaged over if the avg_inconsist argument is set to 1. Otherwise, the simulation will terminate upon detection of such an inconsistency.


----
----
# Sequence File <a name="seq"></a>
The sequence file can specify either the complete sequence using the letters specified in the IDB file, or a repeating pattern. If the full sequence is provided, it must be given in a single line. In case of a repeating pattern, the first line of the sequence file must begin with the word 'repeat', followed by the repeating pattern. For instance, to represent a PolyA sequence, use the following format:
```
repeat
a
```

----
----
# IOPolyMC <a name="iopolymc"></a>
[IOPolyMC](https://github.com/eskoruppa/IOPolyMC) is a python module for reading PolyMC output files in python as well as generating various PolyMC input files. 
Current functionality inclues:
- reading state files (load_state, read_state, read_spec)
- reading and writing xyz files (load_xyz, read_xyz, read_xyz_atomtypes, write_xyz)
- reading and writing IDB files (read_idb,write_idb)
- reading and writing restart files (read_restart, write_restart)
- reading theta files (load_thetas, read_thetas)
- converting state to pdb (state2pdb)
    - This requires base pair discretization with the appropriate intrinsic twist
- reading, writing, and querying input files (read_input, write_input, querysims,simfiles)
- generating a PolyMC configuration by interpolating between given points (pts2config, config2triads, pts2xyz, pts2restart)
    - For now spline interpolation is not included in the package
- finding all unique DNA oligomers of given length N (dna_oligomers)

<!-- ----
----
## Electrostatics
---- -->


<!-- ----
----
## Holonomic constraints
Constraints used in [Ref 4](#vand22).
---- -->

----
# Citing PolyMC <a name="citing"></a>

Please cite [Vanderlinden and Skoruppa et al.](#vand22) for any work that uses the PolyMC simulation package.


----
# Publications <a name="publications"></a>

1. E. Skoruppa, A. Voorspoels, J. Vreede, and E. Carlon. [Length-scale-dependent elasticity in DNA from coarse-grained and all-atom models](https://doi.org/10.1103/PhysRevE.103.042408).
*Phys. Rev. E*, 103:042408, 2021

2. M. Segers, E. Skoruppa, J. A. Stevens, M. Vangilbergen, A. Voorspoels, and E. Carlon. [Comment on “Flexibility of short DNA helices with finite-length effect: From base pairs to tens of base pairs” [J. Chem. Phys. 142, 125103 (2015)]](https://doi.org./10.1063/5.0055349). *J. Chem. Phys.*, 155(2):027101, 2021

3. <a name="skor22"></a>E. Skoruppa and E. Carlon. [Equilibrium fluctuations of DNA plectonemes](https://doi.org/10.1103/PhysRevE.106.024412). *Phys. Rev. E*, 106:024412, 2022

4. <a name="vand22"></a>W. Vanderlinden, E. Skoruppa , P. Kolbeck, E. Carlon, and J. Lipfert. [DNA fluctuations reveal the size and dynamics of topological domains](https://doi.org/10.1093/pnasnexus/pgac268). *PNAS Nexus*, 1:pgac268, 2022

5. <a name="es_phd"></a>E. Skoruppa [Physical Modeling of DNA and DNA-Protein Interactions](https://kuleuven.limo.libis.be/discovery/fulldisplay?docid=lirias3955698&context=SearchWebhook&vid=32KUL_KUL:Lirias&lang=en&search_scope=lirias_profile&adaptor=SearchWebhook&tab=LIRIAS&query=creator%2Cexact%2CU0118787%2CAND&facet=creator%2Cexact%2CU0118787&mode=advanced). PhD thesis, KU Leuven, 2022


----
# References <a name="refs"></a>

6. <a name="klen00"></a>K. Klenin and J. Langowski. [Computation of writhe in modeling of supercoiled DNA](https://doi.org/10.1002/1097-0282(20001015)54:5<307::AID-BIP20>3.0.CO;2-Y). *Biopolymers*, 54(5):307–317, 2000.

7. <a name="full78"></a>F. B. Fuller. [Decomposition of the linking number of a closed ribbon: A problem from molecular biology](https://doi.org/10.1073/pnas.75.8.3557). Proc. Natl. Acad. Sci. U.S.A., 75:3557–3561, 1978.

