# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:02:13 2020
@Description: Convert .csv to ROOT format
@author: benwi
"""
import file_tools as ft

# Make datafiles
ft.convert_all_csv_to_root("/mnt/d/PhD/Work/Majorana/signal_simulation/csv_samples",
                        "/mnt/d/PhD/Work/Majorana/signal_simulation/root_samples",
                        list_files=True)    

ft.make_file_list("../root_samples", "list_of_root_files.txt", ".root", list_files=True)
