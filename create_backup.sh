#!/bin/bash
#Script to backup everything in current directory execpt for directories containing large data files
cwd=${PWD##*/}
dir="../Backups/"$cwd"-"`date +"%d-%m-%Y"`
echo "Making new directory: " $dir
mkdir -p $dir
mkdir -p /mnt/d/PhD/Work/Majorana/$dir
rsync -ah --info=progress2 * ../Backups/$dir --exclude simulation_files --exclude root_samples --exclude hepmc_samples --exclude old_data
rsync -ah --info=progress2 * /mnt/d/PhD/Work/Majorana/Backups/$dir --exclude simulation_files --exclude root_samples --exclude hepmc_samples --exclude old_data
