#!/bin/bash
# Script to automatically run MG from current directory

# Execute MadGraph script
#for script in generator_scripts/*
#do  
#    scriptname=$(basename "${script}")
#    cp $script /home/benwi/MG5_aMC_v2_7_3_py3/ || echo "Command Failed: cp gen_script.txt /home/benwi/MG5_aMC_v2_7_3_py3/"
#    cd /home/benwi/MG5_aMC_v2_7_3_py3 
#    echo "Executing MadGraph with gen_script.txt"
#    time ./bin/mg5_aMC $scriptname || echo "Command Failed: time ./bin/mg5_aMC gen_script.txt"
#done

# Moving files to working directory
cd /mnt/e/MG5_aMC_v2_6_7
echo "Moving files..."
rsync -ah --info=progress2 *NoHiggs /mnt/e/PhD/Work/Majorana/signal_simulation/simulation_files || echo "Command Failed: rsync -ah --info=progress2 *MW* /mnt/d/PhD/Work/Majorana/signal_simulation/simulation_files"
#rm -r *NoHiggs || echo "Command Failed: rm -r *MW*"
cd -

# Decompress hepmc.gz output files and move them to samples directory
for d in simulation_files/*
do  
    echo "d = " $d
    (cd "$d"/Events/run_01 && gunzip tag_1_pythia8_events.hepmc.gz || echo "Command Failed: gunzip tag_1_pythia8_events.hepmc.gz" && rsync -ah --info=progress2 tag_1_pythia8_events.hepmc ../../../../hepmc_samples/$(basename "${d}").hepmc)
done
