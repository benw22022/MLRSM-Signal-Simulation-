# Script to run rivet analysis

# Go to correct directory
cd rivet_analysis

# Build plugin
rivet-buildplugin -j 12 -v MLRSM_mumujj.cc 

# Decompress hepmc.gz output files and move them to samples directory
for d in ../hepmc_samples/*
do  
    fname=$(basename "${d}")
    fnamestrip=${fname%.*}

      echo "d = " $d
      echo "fnamestrip = " $fnamestrip

    rivet --analysis=MLRSM_mumujj --pwd -o whatever.yoda $d &> ../rivet_logs/$fnamestrip.log  && mv rivet_output.txt  ../csv_samples/$fnamestrip.csv
done



