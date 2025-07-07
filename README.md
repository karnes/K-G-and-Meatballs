# K-G and MeatBalls (KGMB)

## Project Description

K-G and MeatBalls (KGMB) was developed to accelerate 
the setup of polymer composite simulations. KGMB accepts 
user input as a text-based input file and generates 
starting conficurations and LAMMPS imput files that perform the 
initial relaxation / equilibration of these configurations. 
As the name suggests, this approach implements the Kremer-Grest 
bead-spring model [Kremer, K.; Grest, G. S. J. Chem. Phys. 1990, 92 (8), 5057–5086. https://doi.org/10.1063/1.458541
] within the LAMMPS molecular dynamics environment. 
KGMB can create configurations with monodisperse, linear polymer 
strands of arbitrary length and the core utility is the ability 
to optionally include spherical filler particles. 

Please cite the work [CITE] when using this code.

## Getting Started

Generation of LAMMPS simulation input files for polymer composites.
*Requirements:* PACKMOL [CITE] installed at a known $PATH
-- default path = ~/packmol/packmol (user can revise in makeLoaded.c)

to use:
1) compile makeLoaded.c
-- gcc -lm makeLoaded.c -o loaded
2) compile makeFibbSphere.c
-- gcc -lm makeFibbSphere.c -o sphere
3) revise runPbatch.template as necessary 
4) revise c shell script makeConfigs.sh to reflect desired inputs
5) run the script 'csh makeConfigs.sh'

## Contributing

K-G and MeatBalls is distributed under the terms of the MIT license. 
All new contributions must be made under this license.

LLNL-CODE-2006739

