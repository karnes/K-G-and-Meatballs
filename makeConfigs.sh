#!/bin/csh
#
foreach beads(100)
foreach rhoStar(87) # /100
foreach packFrac(150) # 250) #350 425 500) # /1000
foreach partRad(2) # 4)
foreach ftemp(0.3) # 1.0)
foreach ljcut(1.22462) # 1.75)
#
mkdir FENE_${beads}chain_${partRad}rad_${packFrac}_T${ftemp}_LJ${ljcut}
cd    FENE_${beads}chain_${partRad}rad_${packFrac}_T${ftemp}_LJ${ljcut}
#
sed -e "s/TEMP/${ftemp}/" ../in.template | sed -e "s/LJ_CUT/${ljcut}/g" > in.template
cp ../sphere .
cp ../makePartBox.inp .
sed -e "s/BEADS_PER_STRAND/${beads}/" ../input.template | sed -e "s/RHO_STAR/${rhoStar}/" | sed -e "s/PACK_FRACTION/${packFrac}/" | sed -e "s/PART_RAD/${partRad}/" > input.rho${rhoStar}
#
../loaded input.rho${rhoStar} 0
#
sed -e "s/JOB/r${rhoStar}eq/" ../runPbatch.template > runPbatch.rho${rhoStar}.equilib
rm slurm*
*sbatch runPbatch.rho${rhoStar}.equilib
#
cd ..
#
end
end
end
end
end
end
