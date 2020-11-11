# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 12:03:45 2020
@ Module containing useful functions for dealing with files
@author: benwi
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:02:13 2020
@Description: Convert .csv to ROOT format
@author: benwi
"""
__name__ = "file_tools"
import uproot
import pandas as pd
import os


def win_path(filepath):
    '''
    Function to force file path to be in windows format

    Parameters
    ----------
    filepath : str
        A file path

    Returns
    -------
    filepath : str
        A file path that should now be in windows format

    '''
    if filepath[:5] == "/mnt/":
        filepath = filepath[5:6] + ":\\" + filepath[7:]         
    filepath = filepath.replace('/', '\\')
    return filepath

def linux_path(filepath):
    '''
    Function to force file path to be in linux format

    Parameters
    ----------
    filepath : str
        a file path

    Returns
    -------
    filepath : str
        A file path that should now be in the linux format

    '''
    
    if filepath[1:3] == ":\\":
        filepath = "/mnt/" + filepath[:1] + filepath[2:]     
    filepath = filepath.replace('\\', '/')
    return filepath
    
def path_convert(filepath, operating_sys_known=False):
    '''
    Function to convert filepaths between linux and windows formats depending
    on os where program is run

    Parameters
    ----------
    filepath : str
        A path to a file/folder
    
    operating_sys (optional): str
        If os is known, no need to perform checks

    Returns
    -------
    filepath : str
        The file path converted to be compatible with current os

    '''
    # If Windows ('nt')
    if os.name == 'nt' and operating_sys_known is False:
        filepath = win_path(filepath)
            
    # If WSL2 Ubuntu ('posix')
    elif os.name == 'posix' and operating_sys_known is False:
        filepath = linux_path(filepath)
    return filepath    
        
    return filepath
    
def csv_to_dataframe(input_file, op_sys_known=None):
    """
    Function to convert a csv file to a pandas dataframe

    Parameters
    ----------
    input_file : str
        Filepath to csv datafile
    op_sys : str, optional
        A. The default is None.

    Returns
    -------
    pandas dataframe
        Returns a pandas dataframe containing the data from the csv
        Returns an empty dataframe if file could not be opened

    """
    data_df = pd.read_csv(path_convert(input_file, operating_sys_known=op_sys_known))

    try:
        data_df = pd.read_csv(path_convert(input_file, operating_sys_known=op_sys_known))
        return data_df
    except:
        print("unable to open file: " + input_file)
        empty_df = pd.DataFrame({'Empty Dataframe' : []})    
        return empty_df


def check_root_filepath(filepath):
    """
    Function to check that outfile is really a path to a .root file

    Parameters
    ----------
    outfile : str
         Path to .root file

    Returns
    -------
    bool
        True if filepath ends in .root, False otherwise

    """
    if ".root" in filepath:
        return True 
    else:
        print("ERROR: " + filepath + " is not a path to a .root file")
        return False
    

def dataframe_to_root_file(data_df, output_file, variables=None,
                           output_tree="output_tree", 
                           additional_data=None):
    """
    Uses upRoot to convert a dataframe into a root file

    Parameters
    ----------
    data_df : pandas dataframe
        Dataframe to be converted into root file
    output_file : str
        Output file name/path
    variables : [str], optional
        List of variables in dataframe to be saved to array. 
        The default is None. If None, all columns of dataframw will be saved
    output_tree : str, optional
        Name of output tree. The default is "output_tree".
    additional_data : pandas dataframe, optional
        An additional dataframe that can be added to root file if required. 
        The default is None

    Returns
    -------
    None.

    """
    # Check that we are actually writing to a .root file, if not return
    if check_root_filepath(output_file) is False:
        return
           
    # Try and save root file
    try:
        with uproot.recreate(output_file, compression=uproot.LZMA(5)) as f:
            
            # Constuct tree template and branch dictionary
            output_arrays = []
            output_branches = [uproot.newbranch("f8")]
            branch_names = []
            tree_name = "output_tree"
                
            # If list of variables is not given- write every column to root file
            if variables is None:
                variables = data_df.columns
            
            # Loop over all columns and copy contents into output arrays 
            for i in range(0, len(variables)):
                output_arrays.append(data_df[variables[i]].to_numpy().ravel().tolist())
                branch_names.append(variables[i])
                output_branches.append(uproot.newbranch("f8"))
            
            # Add any additional data to arrays)
            if isinstance(additional_data, pd.DataFrame):
                add_variables = additional_data.columns
                for i in range(0, len(add_variables)):
                    output_arrays.append(additional_data[add_variables[i]].to_numpy().ravel().tolist())
                    branch_names.append(add_variables[i])
                    output_branches.append(uproot.newbranch("f8"))
                
            # Create dictionaries
            names_and_branches = zip(branch_names, output_branches)
            names_and_data = zip(branch_names, output_arrays)
            branchdict = dict(names_and_branches)
            datadict = dict(names_and_data)
            
            # Create output tree and add branches from dictionary
            tree = uproot.newtree(branchdict, compression=uproot.LZ4(4)) 
            f[tree_name] = tree
            f[tree_name].extend(datadict)
            print("Written: " + output_file)
            
    except:
        print("Unable to create output file: " + output_file)


def csv_to_root_file(input_file, output_file, variables=None,
                      output_tree="output_tree", additional_data=[], op_sys_known=False):
    """
    Function to convert csv file to root file

    Parameters
    ----------
    input_file : str
        DESCRIPTION.
    output_file : str
        Output file name/path
    variables : [str], optional
        List of variables in dataframe to be saved to array. 
        The default is None. If None, all columns of dataframw will be saved
    output_tree : str, optional
        Name of output tree. The default is "output_tree".
    additional_data : pandas dataframe, optional
        An additional dataframe that can be added to root file if required. 
        The default is None
    op_sys_konwn : bool, optional
        Change to True if operating system code is running from is known beforehand.
        If False checks are made to convert between windows and linux filepath convetions. 
        The default is False.

    Returns
    -------
    None.

    """
    # Convert input and output paths if needed
    output_file = path_convert(output_file, op_sys_known)
    input_file = path_convert(input_file, op_sys_known)
    
    if check_root_filepath(output_file) is True:
        data_df = csv_to_dataframe(input_file)
        output_file = path_convert(output_file, op_sys_known)
        
        if data_df.empty is False:
            dataframe_to_root_file(data_df, output_file, variables=None,
                               output_tree="output_tree", additional_data=[])        


def get_files_of_type_in_dir(dir_path, file_type, list_files=False, op_sys=None):
    """
    Function to return a list all files of a certain type in a directory 

    Parameters
    ----------
    dir_path : str
        Directory to search
    file_type : str
        File extension to search for e.g. .txt, .csv, .root etc..
    list_files : Bool, optional
        If true, lists all files found. The default is False.

    Returns
    -------
    file_list : list[str]
        A list of all files of that specific type in that directory
        Returns an empty list if no files found

    """
    dir_path = path_convert(dir_path)
    file_list = []

    for file in os.listdir(dir_path):
        if file.endswith(file_type):
            if op_sys == 'w':
                file_list.append(win_path(os.path.join(dir_path, file)))
            elif op_sys == 'l':
                file_list.append(linux_path(os.path.join(dir_path, file)))
                print(linux_path(os.path.join(dir_path, file)))
            elif op_sys != 'l' or op_sys != 'w':
                file_list.append(os.path.join(dir_path, file))
    if list_files is True:
        print("Found: " + str(len(file_list)) + " files")
        for file in os.listdir(dir_path):
            if file.endswith(file_type):
                print(os.path.join(dir_path, file))
    
    if not file_list:
        print("No " + file_type + " were found in " + dir_path)
        
    return file_list


def convert_all_csv_to_root(input_dir_path, output_dir_path, list_files=False):
    """
    Function to convert all .csv files in a directory to root files
    Default root file configureation used

    Parameters
    ----------
    input_dir_path : str
        Directory to file containing the csv files to convert
    output_dir_path : str
        Path to destination directory

    Returns
    -------
    None.

    """
    file_list = get_files_of_type_in_dir(input_dir_path, ".csv", list_files)
    for file in file_list:
        output_file = path_convert(output_dir_path)
        if os.name == 'nt':
            index = file.rfind('\\')
            output_file = output_file + file[index:-4] + ".root"    
        else:    
            index = file.rfind('/')
            output_file = output_file + file[index:-4] + ".root"    
        
        csv_to_root_file(file, output_file)


def write_list_to_file(outfile_name, data_list, delimiter="\n"):
    '''
    Function to take a list and write it to a file. Elements seperated by a
    particular delimiter

    Parameters
    ----------
    outfile_name : str
        The name/path to the output file to be written
    data_list : []
        A list of data to be written to file
    delimiter : char, optional
        Character to delimit element in list. The default is "\n".

    Returns
    -------
    None.

    '''
    file = open(outfile_name, "w")
    for element in data_list:
        file.write(str(element) + delimiter)


def make_file_list(dir_path, outfile_name, file_type, list_files=False, op_sys='l'):
    '''
    Function to make a file listing all files of a particular type in a directory

    Parameters
    ----------
    dir_path : str
        File path to directory which contains the files to be listed
    outfile_name : str
        Name/path to file where list of files will be written
    file_type : str
        File extension to search for 
    list_file : Bool, optional
        List the files found in directory in the terminal. The default is False.
    op_sys : char, optional
        Operating system format for filepath to be written to. 
        The default is 'L' for linux. Other option is 'W' for windows

    Returns
    -------
    None.

    '''
    dir_path = path_convert(dir_path)
    outfile_name = path_convert(outfile_name)
    files = get_files_of_type_in_dir(dir_path, file_type, list_files, op_sys='l')
    print("\n   **** \n ")
    print(files)
    write_list_to_file(outfile_name, files)