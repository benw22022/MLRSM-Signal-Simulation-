# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 12:05:44 2020
@Script to automatically produce generation script files
@author: benwi
"""

import pandas as pd
import math
import file_tools as ft


class Generator_Script:
        
    def __init__(self, region, top_lines, filename):
        self.m_region = region
        self.m_lines = top_lines
        self.m_filename = filename
        self.m_number_of_files = 0
        self.m_output_files = []
                
    def __template_lines(self, set_MW2, set_MN5, set_MN4=10000, set_MN6=10000,
                       set_WW2="auto", set_WN4="auto", set_WN5="auto", set_WN6="auto", calc_widths=True):
        
        if calc_widths is True:
            lines = [self.m_output_files[len(self.m_output_files)-1],
                    "launch", 
                    "shower=pythia8", 
                    "set MW2" + " = " + str(set_MW2), 
                    "set MN4" + " = " + str(set_MN4), 
                    "set MN5" + " = " + str(set_MN5),
                    "set MN6" + " = " + str(set_MN6),
                    "set WW2" + " = " + str(set_WW2),
                    "set WN4" + " = " + str(set_WN4),
                    "set WN5" + " = " + str(set_WN5),
                    "set WN6" + " = " + str(set_WN6)]
            return lines
        lines = [self.m_output_files[len(self.m_output_files)-1],
                    "launch", 
                    "shower=pythia8", 
                    "set MW2" + " = " + str(set_MW2), 
                    "set MN4" + " = " + str(set_MN4), 
                    "set MN5" + " = " + str(set_MN5),
                    "set MN6" + " = " + str(set_MN6)]
        return lines


    def add_mass_point(self, region, WR_mass, N5_mass,
                     process="pp_to_mumujj", N46_masses="10", mode="WRWR", higgs_opt="NoHiggs",
                     calc_widths=True):
        
        self.m_output_files.append("output " + region + "_" + process +"_MWR_" + str(WR_mass) + \
            "_MN5_" + str(N5_mass) + "_MN46_" + N46_masses + "_" + mode + "_" + higgs_opt + \
                "_CalcWidths_" + str(calc_widths))
            
        self.m_lines = self.m_lines + self.__template_lines(WR_mass*1000, N5_mass*1000, calc_widths=True)
        self.m_number_of_files += 1
    
    
    def add_lines_to_file(self, csv_file, limit=0): 
        Mos_region = ft.csv_to_dataframe(ft.path_convert(csv_file))
        models = Mos_region["Model"].to_numpy()
        model_count = 0
        for model in models:
            if "WRWR" in model:
                print("MWR = " + model[3:6] + "  MN5 = " + model[10:13])
                self.add_mass_point(self.m_region, WR_mass=float(model[3:6]), 
                                             N5_mass = float(model[10:13]))
                model_count += 1
                if limit != 0 and model_count >= limit:
                    break

    def write_to_file(self, verbose=True):
        ft.write_list_to_file(self.m_filename, self.m_lines)
        print("Written %i mass points to %s" % (self.m_number_of_files, self.m_filename))


# Main code body
if __name__ == "__main__":
    
    # Lines at the top of every script
    top_lines = ["import model MLRSM",
                 "define bsmhiggs h h2 h02 h03 h+ hp2 hl++ hr++ h3 a02 h- hm2 hl-- hr--",
                 "generate p p > mu+ mu+ j j / bsmhiggs w+ w- @1",
                 "add process p p > mu- mu- j j / bsmhiggs w+ w- @2"]
    
    
    # Region A1
    A1_gen_script = Generator_Script("A1", top_lines, ft.path_convert("../generator_scripts/A1_gen_script.txt"))
    A1_gen_script.add_lines_to_file("/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/A1_xs_allbkgs.csv", limit=3)
    A1_gen_script.write_to_file()
    
    # Region A2
    A2_gen_script = Generator_Script("A2", top_lines, ft.path_convert("../generator_scripts/A2_gen_script.txt"))
    A2_gen_script.add_lines_to_file("/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/A2_xs_allbkgs.csv", limit=3)
    A2_gen_script.write_to_file()

    # Region B1
    B1_gen_script = Generator_Script("B1", top_lines, ft.path_convert("../generator_scripts/B1_gen_script.txt"))
    B1_gen_script.add_lines_to_file("/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/B1_xs_allbkgs.csv", limit=3)
    B1_gen_script.write_to_file()
    
    # Region B2
    B2_gen_script = Generator_Script("B2", top_lines, ft.path_convert("../generator_scripts/B2_gen_script.txt"))
    B2_gen_script.add_lines_to_file("/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/B2_xs_allbkgs.csv", limit=3)
    B2_gen_script.write_to_file()

    # Region D
    D_gen_script = Generator_Script("D", top_lines, ft.path_convert("../generator_scripts/D_gen_script.txt"))
    D_gen_script.add_lines_to_file("/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/D_xs_allbkgs.csv", limit=3)
    D_gen_script.write_to_file()
    
    print("#######################################################")
    print("Total number of files: ")
    total_number_of_files = (A1_gen_script.m_number_of_files +
                            A2_gen_script.m_number_of_files +
                            B1_gen_script.m_number_of_files +
                            B2_gen_script.m_number_of_files +
                            D_gen_script.m_number_of_files)
    total_time_hrs = total_number_of_files * 3
    total_time_days_hrs = [math.floor(total_time_hrs / 24), total_time_hrs % 24]
    print("Region A1: %i" % A1_gen_script.m_number_of_files)
    print("Region A2: %i" % A2_gen_script.m_number_of_files)
    print("Region B1: %i" % B1_gen_script.m_number_of_files)
    print("Region B2: %i" % B2_gen_script.m_number_of_files)
    print("Region D: %i" % D_gen_script.m_number_of_files)
    print("\nTotal Number of files = %i" % (total_number_of_files))
    print("Assuming ~3hrs per file this will take: %i days %i hrs" %(total_time_days_hrs[0],
                                                                     total_time_days_hrs[1]))
    print("#######################################################")





































'''
# Region A1
WR_masses = [1000, 1500, 2000, 2500, 3000, 3500, 4000, 6000, 7000]
N5_masses = [[2000], [2000, 2500, 3000], [2000, 2500, 3000, 3500], [2500, 3000, 3500], [3000, 3500, 4000], 
             [3500, 4000], [4000], [6000], [6000]]    
gen_script_A1_lines = region_lines("A1", WR_masses, N5_masses, top_lines)
ft.write_list_to_file("../generator_scripts/A1_gen_script.txt", gen_script_A1_lines[0])

# Region A2
WR_masses = [2000, 2000, 3000, 4000]
N5_masses = [[5000] ,[6000] ,[6000], [7000]]
gen_script_A2_lines = region_lines("A2", WR_masses, N5_masses, top_lines)
ft.write_list_to_file("../generator_scripts/A2_gen_script.txt", gen_script_A2_lines[0])

# Region B1
WR_masses = [3500, 4000, 7000]
N5_masses = [[2500], [3000], [5500]]
gen_script_B1_lines = region_lines("B1", WR_masses, N5_masses, top_lines)
ft.write_list_to_file("../generator_scripts/B1_gen_script.txt", gen_script_B1_lines[0])

# Region B2
WR_masses = [5000, 6000, 7000]
N5_masses = [[1500], [1500, 2000], [3500]]
gen_script_B2_lines = region_lines("B2", WR_masses, N5_masses, top_lines)
ft.write_list_to_file("../generator_scripts/B2_gen_script.txt", gen_script_B2_lines[0])

# Region ATLAS
WR_masses = [1000, 1500, 2000, 2500, 3000, 3500, 4000]
N5_masses = [[1000], [1000, 1500], [1000, 1500], [1500, 2000], [1000, 1500, 2000, 2500], [1500, 2000, 2500],
             [1000, 2000, 3000]] 
gen_script_ATLAS_lines = region_lines("ATLAS", WR_masses, N5_masses, top_lines)
ft.write_list_to_file("../generator_scripts/ATLAS_gen_script.txt", gen_script_ATLAS_lines[0])


regions = ["A1", "A2", "B1", "B2", "ATLAS"]
files = [gen_script_A1_lines[1], gen_script_A2_lines[1], gen_script_B1_lines[1], gen_script_B2_lines[1],
         gen_script_ATLAS_lines[1]]
print("\n ***********************************************")
for i in range(0, len(regions)):
    print("Number of files in " + regions[i] +" = %i" % files[i])
print("Number of files in total = " + str(sum(files)))
print("Assuming 15 mins per file, this will take ~%6.2f hours" %(sum(files)/4))
'''

