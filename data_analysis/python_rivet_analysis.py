# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 12:50:35 2020

@author: benwi
"""

from multiprocessing import Pool
import subprocess
import file_tools as ft

def run_rivet_analysis(files):
    
    infile = files[0]
    outfile = files[1]
    
    analysis_result = subprocess.run(
        ["rivet", "--analysis=HeavyN_pp_to_mumujj", "--pwd", "-o", "whatever.yoda", infile 
        ,"&>", "../rivet_logs/$fnamestrip.log", "&& mv rivet_output.txt", outfile]
        , stdout=subprocess.PIPE)
    return analysis_result
    
if __name__ == "__main__":
    
    """
    os.chdir("rivet_analysis")
    compile_result = subprocess.run(["rivet-buildplugin", "-j", "12", "-v", 
                                     "HeavyN_pp_to_mumujj.cc"], stdout=subprocess.PIPE)
    
    print("*****************************************************")
    print ("Compile Return code: %i" % compile_result.returncode)
    print ("Compile: %s" % compile_result.stdout)
    print("*****************************************************\n")
    """
    hepmc_file_list = ft.get_files_of_type_in_dir("../hepmc_samples", ".hepmc")
    
    csv_file_list = []
    for file in hepmc_file_list:
        fname_split = file.split(ft.path_convert("/"))
        csv_file_list.append(ft.path_convert("../csv_samples/") + fname_split[2].strip(".hepmc") + ".csv")
    
    files = [hepmc_file_list, csv_file_list]
    
    agents = 12
    chunksize = 3
    
    print("Analysing %i with %i agents, chunksize = %i" %(len(hepmc_file_list), agents, chunksize))
    
    with Pool(processes=agents) as pool:
        result = pool.map(run_rivet_analysis, files, chunksize)

    # Output the result
    print ('Result:  ' + str(result))
    
    
