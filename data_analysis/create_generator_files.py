# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 12:05:44 2020
@Script to automatically produce generation script files
@author: benwi
"""

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
        """
        Template lines for adding a new mass point 

        Parameters
        ----------
        set_MW2 : float
            mass of WR in GeV
        set_MN5 : float
            mass of N5 in GeV
        set_MN4 : float, optional
            mass of N4 in GeV. The default is 10000.
        set_MN6 : float, optional
            mass of N6 in GeV. The default is 10000.
        set_WW2 : float/str, optional
            Width of WR. The default is "auto".
        set_WN4 : float/str, optional
            Width of N4. The default is "auto".
        set_WN5 : float/str, optional
            Width of N5. The default is "auto".
        set_WN6 : float/str, optional
            Width of N6. The default is "auto".
        calc_widths : bool, optional
            calulate widths. The default is True.

        Returns
        -------
        lines : [str]
            Array of lines to add to file.

        """
        
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
                     process="pp_to_mumujj", N46_masses=10, mode="WRWR",
                     codes={"CalcWidth":1, "HiggsExcl":1, "JetDef":1}):
        '''
        Adds a mass point

        Parameters
        ----------
        region : str
            Lables mass region
        WR_mass : float
            WR mass in GeV
        N5_mass : float
            N5 mass in GeV
        process : str, optional
            lables production process. The default is "pp_to_mumujj".
        N46_masses : float, optional
            masses of N4 and N6. The default is 10.
        mode : float, optional
            production mode. The default is "WRWR".
        higgs_opt : str, optional
            lables Higgs production parameters. The default is "NoHiggs".
        calc_widths : bool, optional
            run width calculation. The default is True.

        Returns
        -------
        None.

        '''
        
        run_codes = "W"+str(codes["CalcWidth"])+"H"+str(codes["HiggsExcl"])+\
                    "J"+str(codes["JetDef"])
        
        self.m_output_files.append("output " +"simulation_files/"+ region + "_" + process +"_MWR_" + str(WR_mass) + \
            "_MN5_" + str(N5_mass) + "_MN46_" + str(N46_masses) + "_" + mode+"_"+run_codes)
            
        self.m_lines = self.m_lines + self.__template_lines(WR_mass*1000, N5_mass*1000, 
                                                            set_MN4=N46_masses*1000, 
                                                            set_MN6=N46_masses*1000, 
                                                            calc_widths=True)
        self.m_number_of_files += 1
    
    
    def add_lines_to_file(self, mass_array=None, csv_file=None, MN46_same=False, N46_mass=10, limit=0): 
        
        if csv_file is not None:
            Mos_region = ft.csv_to_dataframe(ft.path_convert(csv_file))
            models = Mos_region["Model"].to_numpy()
            model_count = 0
            for model in models:
                if "WRWR" in model:
                    print("MWR = " + model[3:6] + "  MN5 = " + model[10:13])
                    
                    if MN46_same == True:
                        N46_mass = float(model[10:13])
                    
                    self.add_mass_point(self.m_region, WR_mass=float(model[3:6]), 
                                                 N5_mass = float(model[10:13]),
                                                 N46_masses = N46_mass)
                    model_count += 1
                    if limit != 0 and model_count >= limit:
                        break
                    
        elif mass_array is not None:
            if MN46_same == True:
                for mass_point in mass_array:
                    self.add_mass_point(self.m_region, WR_mass=mass_point[0], 
                                        N5_mass = mass_point[1],
                                        N46_masses = mass_point[1])
            else:
                for mass_point in mass_array:
                    self.add_mass_point(self.m_region, WR_mass=mass_point[0], 
                                        N5_mass = mass_point[1],
                                        N46_masses = N46_mass)
            

    def write_to_file(self, verbose=True):
        '''
        Writes the lines to file

        Parameters
        ----------
        verbose : bool, optional
            prints info if true. The default is True.

        Returns
        -------
        None.

        '''
        ft.write_list_to_file(self.m_filename, self.m_lines)
        print("Written %i mass points to %s" % (self.m_number_of_files, self.m_filename))


# Main code body
if __name__ == "__main__":
    
    # Lines at the top of every script
    top_lines = ["import model MLRSM",
                 "define bsmhiggs h h2 h02 h03 h+ hp2 hl++ hr++ h3 a02 h- hm2 hl-- hr--",
                 "define j = g u c d s b t u~ c~ d~ s~ b~ t~",
                 "generate p p > mu+ mu+ j j / bsmhiggs w+ w- @1",
                 "add process p p > mu- mu- j j / bsmhiggs w+ w- @2"]
    
    
    # Region A1
    A1_gen_script = Generator_Script("A1", top_lines, ft.path_convert("../generator_scripts/A1_gen_script.txt"))
    A1_gen_script.add_lines_to_file(csv_file="/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/A1_xs_allbkgs.csv", limit=3, MN46_same=True)
    A1_gen_script.write_to_file()
    
    # Region A2
    A2_gen_script = Generator_Script("A2", top_lines, ft.path_convert("../generator_scripts/A2_gen_script.txt"))
    A2_gen_script.add_lines_to_file(csv_file="/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/A2_xs_allbkgs.csv", limit=3, MN46_same=True)
    A2_gen_script.write_to_file()

    # Region B1
    B1_gen_script = Generator_Script("B1", top_lines, ft.path_convert("../generator_scripts/B1_gen_script.txt"))
    B1_gen_script.add_lines_to_file(csv_file="/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/B1_xs_allbkgs.csv", limit=3, MN46_same=True)
    B1_gen_script.write_to_file()
    
    # Region B2
    B2_gen_script = Generator_Script("B2", top_lines, ft.path_convert("../generator_scripts/B2_gen_script.txt"))
    B2_gen_script.add_lines_to_file(csv_file="/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/B2_xs_allbkgs.csv", limit=3, MN46_same=True)
    B2_gen_script.write_to_file()

    # Region D
    D_gen_script = Generator_Script("D", top_lines, ft.path_convert("../generator_scripts/D_gen_script.txt"))
    D_gen_script.add_lines_to_file(csv_file="/mnt/d/PhD/Work/Majorana/Mos_MPhys/Data/D_xs_allbkgs.csv", limit=3, MN46_same=True)
    D_gen_script.write_to_file()
    
    # Region Test
    test_mass_points = [[1,1],
                        [1,2],
                        [2,1],
                        [2,2],
                        [2,3],
                        [3,2],
                        [3,3]]
    
    Test_gen_script = Generator_Script("Test", top_lines, ft.path_convert("../generator_scripts/test_gen_script.txt"))
    Test_gen_script.add_lines_to_file(mass_array=test_mass_points, MN46_same=True)
    Test_gen_script.write_to_file()
    
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
    Codes:
        Higgs excluded- 0:No  1:Yes
        Widths calculated- 0:No 1:Yes
        Jet definition- 0:MG default  1:Inc bbar, ttbar
    
    '''