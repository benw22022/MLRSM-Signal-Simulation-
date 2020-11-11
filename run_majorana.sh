 # Script to run rivet analysis

 # Go to correct directory
 cd ~/mnt/
 source go_to_workDir.sh
 cd signal_simulation/rivet_analysis

 # Run aliases
 -v $PWD:$PWD -w $PWD -u `id -u $USER`:`id -g`
 alias rivet='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet'
alias rivet-mkanalysis='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet-mkanalysis'
alias rivet-buildplugin='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet-buildplugin'
alias rivet-mkhtml='docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet-mkhtml'
alias yodamerge='docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 yodamerge'

# Build plugin
rivet-buildplugin majorana_test.cc

# Run analysis
cd ../hepmc_samples
rivet --analysis=majorana_test --pwd -o tag_1_pythia8_events.hepmc &> ../rivet_logs/MW_1_MN_2_WRWR.log  && mv <txt_file>.txt  ../csv_samples/MW_1_MN_2_WRWR.csv 