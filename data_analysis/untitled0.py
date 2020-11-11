# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:53:08 2020
@description: create easy to read .csv files from rivet logs
@author: benwi
"""
import file_tools as ft
import pandas as pd

class Entry:
    """
    Class to hold log entry data
    """
    def __init__(self, wr_mass, n5_mass, mass_region):
        self.WR_mass = wr_mass
        self.N5_mass = n5_mass
        self.region = mass_region
        self.rivet_xsec = -999
        self.rivet_xsec_err = -999
        self.rivet_sum_of_weights = -999
        self.rivet_sum_of_weights_err = -999
        self.mg_xsec = -999
        self.mg_xsec_err = -999
        self.mg_sum_of_weights = -999
        self.mg_sum_of_weights_err = -999
        self.extra_info = ""

    def set_rivet_data(self, xsec, xsec_err, sow, sow_err):
        self.rivet_xsec = xsec
        self.rivet_xsec_err = xsec_err
        self.rivet_sum_of_weights = sow
        self.rivet_sum_of_weights_err = sow_err

    def set_mg_data(self, xsec, xsec_err, sow, sow_err):
        self.mg_xsec = xsec
        self.mg_xsec_err = xsec_err
        self.mg_sum_of_weights = sow
        self.mg_sum_of_weights_err = sow_err

    def set_extra_info(self, extra_info):
        self.extra_info = extra_info


class Entries:
    """
    Class to hold a list of Entry classes
    Methods allow it to return a list of data of a particular atrribute
    In hindsight probably should of just used a dataframe lol ;)
    """

    def __init__(self):
        self.entries = []

    def append(self, new_entry):
        if type(new_entry) is Entry:
            self.entries.append(new_entry)

    def array(self, var):
        arg_arr = []
        for entry_in_list in self.entries:
            arg_arr.append(vars(entry_in_list).get(var))
        return arg_arr



if __name__ == "__main__":

    log_file_path = ft.path_convert("../rivet_logs")
    log_file_list = ft.get_files_of_type_in_dir(log_file_path, ".log")

    print("Found: %i log files" % len(log_file_list))

    entries = Entries()

    for filename in log_file_list:

        WR_mass = filename[int(filename.find("MWR")+4): int(filename.find("MWR")+7)]
        N5_mass = filename[int(filename.find("MN5")+4): int(filename.find("MN5")+7)]

        region = ""
        if filename[:1] != "D":
            region = filename[int(filename.find("_pp")-2): int(filename.find("_pp"))]
        else:
            region = "D"

        entry = Entry(WR_mass, N5_mass, region)

        read_file = False
        try:
            float(WR_mass)
            float(N5_mass)
            read_file = True
        except ValueError:
            print("Could not access mass data for" + filename)

        if read_file:
            with open(filename, 'r') as file:
                data = file.readlines()
                for line in data:
                    if "Total XS MadGraph:" in line:
                        entry.mg_xsec = line.replace("Total XS MadGraph:", '').replace(" ", '').replace("\n", '')
                    elif "Total XS error:" in line:
                        entry.mg_xsec_err = line.replace("Total XS error:", '').replace(" ", '').replace("\n", '')
                    elif "Total Sum of Weights:" in line and entry.mg_sum_of_weights == -999:
                        entry.mg_sum_of_weights = line.replace("Total Sum of Weights:", '').replace(" ", '').replace("\n", '')
                    elif "Error on Total Sum of Weights:" in line:
                        entry.mg_sum_of_weights_err = line.replace("Error on Total Sum of Weights:", '').replace(" ", '').replace("\n", '')

                    elif "Fiducial crossection:" in line:
                        entry.rivet_xsec = line.replace("Fiducial crossection:", '').replace(" ", '').replace("\n", '')
                    elif "Error of fiducial crossection:" in line:
                        entry.rivet_xsec_err = line.replace("Error of fiducial crossection:", '').replace(" ", '').replace("\n", '')
                    elif "Selected Events Sum of Weights:" in line and entry.rivet_sum_of_weights == -999:
                        entry.rivet_sum_of_weights = line.replace("Selected Events Sum of Weights:", '').replace(" ", '').replace("\n", '')
                    elif "Error on Selected Events Sum of Weights:" in line:
                        entry.rivet_sum_of_weights_err = line.replace("Error on Selected Events Sum of Weights:", '').replace(" ", '').replace("\n", '')

            entry.set_extra_info(filename[filename.find("MN5")+8:].replace(".log", ""))
            entries.append(entry)

    meta_data_dict = {"region" : entries.array("region"),
                      "MWR" : entries.array("WR_mass"),
                      "MN5" : entries.array("N5_mass"),
                      "rivet xsec" : entries.array("rivet_xsec"),
                      "rivet sum of weights" : entries.array("rivet_sum_of_weights"),
                      "mg xsec" : entries.array("mg_xsec"),
                      "mg sum of weights" : entries.array("mg_sum_of_weights"),
                      "rivet xsec err" : entries.array("rivet_xsec_err"),
                      "rivet sum of weights err" : entries.array("rivet_sum_of_weights_err"),
                      "mg xsec err" : entries.array("mg_xsec_err"),
                      "mg sum of weights err" : entries.array("mg_sum_of_weights_err"),
                      "extra info" : entries.array("extra_info")}

    meta_data_df = pd.DataFrame.from_dict(meta_data_dict)

    try:
        meta_data_df.to_csv(ft.path_convert("../rivet_logs/log_data.csv"), index=False)
    except PermissionError:
        print("File is already open! Close file and try again")
