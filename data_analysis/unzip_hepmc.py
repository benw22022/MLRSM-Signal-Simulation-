# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 11:40:54 2020
@Unzip all hepmc.gz files and place them in a folder
@Jobs will be executed in parallel
@author: benwi
"""

from multiprocessing import Pool
from time import sleep
import file_tools as ft
import gzip 
import shutil
import os

def decompress_file(files):
    '''
    Function to decompress .gz files

    Parameters
    ----------
    files : [str]
        Two-componant array- first is a list of input file paths,
        second is a list of output paths. Done this way to make the multiproccessing
        easier

    Returns
    -------
    None.

    '''
    infile = files[0]
    outfile = files [1]
    
    print(infile)
    
    try:
        with gzip.open(infile, 'rb') as f_in:
            with open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print("Decompressed " + infile)
    except FileNotFoundError as file_error:
        print("File " + infile + " not found")
    
if __name__ == "__main__":
     
    path = ft.path_convert("../simulation_files")    
    infiles = [f.path for f in os.scandir(path) if f.is_dir()]
    outfiles = []
    for file in infiles:
        fname_split = file.split(ft.path_convert("/"))
        outfiles.append(ft.path_convert("../hepmc_samples/") + fname_split[2] + ".hepmc")

    for i in range(0, len(infiles)):
        infiles[i] = infiles[i] + ft.path_convert("/Events/run_01/tag_1_pythia8_events.hepmc.gz")
    
    files = [infiles, outfiles]
    
    agents = 12
    chunksize = 3
    
    print("Decompressing %i with %i agents, chunksize = %i" %(len(infiles), agents, chunksize))
    
    with Pool(processes=agents) as pool:
        result = pool.map(decompress_file, files, chunksize)

    # Output the result
    print ('Result:  ' + str(result))
