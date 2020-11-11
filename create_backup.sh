#!/bin/bash
#Script to backup everything in current directory execpt for directories containing large data files
cwd=${PWD##*/}
dir="../Backups/"$cwd"-"`date +"%d-%m-%Y"`
echo "Making new directory: " $dir
mkdir -p $dir
rsync -ah --info=progress2 * ../Backups/$dir --exclude simulation_files --exclude root_samples --exclude hepmc_samples --exclude old_data
