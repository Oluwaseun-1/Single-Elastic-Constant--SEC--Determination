#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 15:25:41 2025

@author: marowolo
"""
import os
import sys
import time
import re
import numpy as np
import matplotlib.pyplot as plt
#import pygad
from scipy import stats
import pandas as pd
from numpy.linalg import inv
import random

start = time.time()



file_path =  "E:/New_PhD_Dell_19-03-2025/Result/New Petra/Ti6242-4mm-840C-103S-V2/2D-Integration/"

stress_strain_file = file_path + "Ti64-3mm-750C-103S-V2.xlsx"

micro_strain_file = file_path + ""

strain_book = 'Dspacing_FIT_Elastic' #'Dspacing_FIT'
stress_book = 'Stress-exp-elastic' #'Stress-exp'


sheet_name = 'Dspacing'

a_1 = ["Alpha -{10-10}", "Alpha -{10-10}_90"]
a_2 = ["Alpha -{0002}", "Alpha -{0002}_90"]
a_3 = ["Alpha -{10-11}", "Alpha – {10-11}_90"]
a_4 = ["Alpha -{10-12}", "Alpha – {10-12}_90"]
a_5 = ["Alpha -{11-20}", "Alpha -{11-20}_90"]

b_1 = ["Beta -{110}", "Beta -{110}_90"]
b_2 = ["Beta -{200}", "Beta -{200}_90"]


alpha_planes_dic = {'a_1': a_1, 'a_2': a_2, 'a_3': a_3,
                    'a_4': a_4, 'a_5': a_5}
beta_planes_dic  = {'b_1': b_1, 'b_2': b_2}

alpha_beta_planes_dic = {'a_1': a_1, 'a_2': a_2, 'b_1': b_1, 'a_3': a_3,
                    'a_4': a_4, 'a_5': a_5, 'b_2': b_2}


def get_file(path):
    
    exc_file_name = []

    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.xlsx'):
                exc_file_name.append(os.path.join(root, file))
    
    deg_0_file = deg_90_file = deg_180_file = deg_270_file = None
    
    pattern = re.compile(r'(\d+)D', re.IGNORECASE)
    
    for excel_file in exc_file_name:
        match = re.search(pattern, os.path.basename(excel_file))
        if match:
            id_ = int(match.group(1))  # Extract the numeric angle

            if id_ == 0:
                deg_0_file = excel_file
            elif id_ == 90:
                deg_90_file = excel_file
            elif id_ == 180:
                deg_180_file = excel_file
            elif id_ == 270:
                deg_270_file = excel_file
            else:
                print(f'{excel_file} has no known direction')
        else:
            print(f'{excel_file} does not match the expected pattern')

    return deg_0_file, deg_90_file, deg_180_file, deg_270_file
    

def get_direc_file_head(file_inp, sheet_name):
    
    file_read = pd.read_excel(file_inp, sheet_name=sheet_name.strip(), engine='openpyxl').replace(['', '-', 'NA'], np.nan).astype(float)
    
    header_list = file_read.head()
    
    pattern = re.compile(r'^A|B')
    
    alpha_beta_data = [head for head in header_list if re.search(pattern, head)]
    
    return alpha_beta_data

def oper_direc_file(path, sheet_name, alpha_beta_planes_dic):
    # Get the directional files
    deg_0_file, deg_90_file, deg_180_file, deg_270_file = get_file(path)
    
    deg_0_file_red,  deg_90_file_red  = None, None
    
    if any(file is None for file in [deg_0_file, deg_90_file, deg_180_file, deg_270_file]):
        print('One of deg_0_file, deg_90_file, deg_180_file, deg_270_file is None')
        sys.exit()

    count = 0

    for file_deg in [deg_0_file, deg_90_file, deg_180_file, deg_270_file]:
        if file_deg is None:
            print(file_deg)
            print(f'{file_deg} is empty')
            continue  


        alpha_beta_data = get_direc_file_head(file_deg, sheet_name)

        file_deg = pd.read_excel(file_deg, sheet_name=sheet_name.strip(), engine='openpyxl').replace(['', '-', 'NA'], np.nan).astype(float)

        file_deg = pd.DataFrame(file_deg)

        pattern_id_index_0, pattern_id_index_180 = [],[]
        pattern_id_index_90, pattern_id_index_270 = [], [] 
        
        count += 1

        if all(header in  list(alpha_beta_planes_dic.get(i)[0] for i in alpha_beta_planes_dic.keys()) for header in alpha_beta_data):

            if count == 1:
                deg_0_file = file_deg

            elif count == 3:
                deg_180_file = file_deg

        elif all(header in  list(alpha_beta_planes_dic.get(i)[1] for i in alpha_beta_planes_dic.keys()) for header in alpha_beta_data):

            if count == 2:
                deg_90_file = file_deg
            elif count == 4:
                deg_270_file = file_deg
        else:
            print(f'{os.path.basename(file_deg)} is unknown')

    if len(deg_0_file) == len(deg_180_file) and deg_0_file['Pattern-id'].all() == deg_180_file['Pattern-id'].all():

        
        avg_ = {alpha_beta_planes_dic[key][0]: [] for key in alpha_beta_planes_dic}
        
        count = 0

        for header in avg_.keys():
            count +=1
            deg_0 = deg_0_file[header] 
            deg_180 = deg_180_file[header]

  
            for data_0_val, data_180_val in zip(deg_0, deg_180):

                if 0 in (data_0_val, data_180_val):

                    avg_[header].append(data_0_val if data_0_val != 0 else data_180_val if data_180_val != 0 else 0)

                else:

                    avg_[header].append(0.5*(data_0_val + data_180_val))

    else:
        
        for (index_0, index_180), (pattern_id_0, pattern_id_180) in zip(enumerate(deg_0_file['Pattern-id']), enumerate(deg_0_file['Pattern-id'])):
            
            if pattern_id_0 != pattern_id_180:
                
                pattern_id_index_0.append(index_0)
                pattern_id_index_180.append(index_180)
        pattern_id_index_0_180 = list(set(pattern_id_index_0).union(pattern_id_index_180))
        
        deg_0_file_red = deg_0_file.pop(pattern_id_index_0_180)
        
        deg_180_file_red = deg_180_file.pop(pattern_id_index_0_180)
        
        for header in len(alpha_beta_planes_dic):
            deg_0 = deg_0_file_red[header] 
            deg_180 = deg_180_file_red[header]
            
            for data_0_val, data_180_val in zip(deg_0, deg_180):

                if 0 in (data_0_val, data_180_val):

                    avg_[header].append(data_0_val if data_0_val != 0 else data_180_val if data_180_val != 0 else 0)

                else:

                    avg_[header].append(0.5*(data_0_val + data_180_val))
                    

    if len(deg_90_file) == len(deg_270_file) and deg_90_file['Pattern-id'].all() == deg_270_file['Pattern-id'].all():

        avg_90_ = {alpha_beta_planes_dic[key][1]: [] for key in alpha_beta_planes_dic}

        for header_90 in avg_90_.keys():
            #counter += 1
           
            deg_90 = deg_90_file[header_90] 
            deg_270 = deg_270_file[header_90]
            
            for data_90_val, data_270_val in zip(deg_90, deg_270):
                    
                    if 0 in (data_90_val, data_270_val):

                        avg_90_[header_90].append(data_90_val if data_90_val != 0 else data_270_val if data_270_val != 0 else 0)
                        
                    else:

                        avg_90_[header_90].append(0.5*(data_90_val + data_270_val))
                     
    else:
        
        for (index_90, index_270) , (pattern_id_90, pattern_id_270) in zip(enumerate(deg_90_file['Pattern-id']), enumerate(deg_270_file['Pattern-id'])):
            
            if pattern_id_90 != pattern_id_270:
                
                pattern_id_index_90.append(index_90)
                pattern_id_index_270.append(index_270)
                
                print('pattern_id_index_90, pattern_id_index_270', pattern_id_index_90, pattern_id_index_270)
                
        pattern_id_index_90_270 = list(set(pattern_id_index_90).union(pattern_id_index_270))
        
        deg_90_file_red = deg_90_file.pop(pattern_id_index_90_270)
        
        deg_270_file_red = deg_270_file.pop(pattern_id_index_90_270)
        
        for header_90 in alpha_beta_planes_dic:
            deg_90 = deg_90_file_red[header_90] 
            deg_270 = deg_270_file_red[header_90]
            
            for (index_90, index_270), (data_90, data_270) in zip(enumerate(deg_90), enumerate(deg_270)):
                
                if 0 in (data_90, data_270):
                    
                    avg_90_[header_90].append(data_90_val if data_90_val != 0 else data_270_val if data_270_val != 0 else 0)
                else:
                    avg_90_[header_90].append(0.5*(data_90_val + data_270_val))
    
    if deg_0_file_red and deg_90_file_red is not None:
                    
        new_0_180_avg = {'Pattern-id': deg_0_file_red['Pattern-id'], 
                         'Macro-True Strain':deg_0_file_red['Macro-True Strain'], 
                          **{header: avg_[header] for header in avg_} }


        new_90_270_avg = {'Pattern-id':deg_90_file_red['Pattern-id'],
                        'Macro-True Strain_90':deg_90_file_red['Macro-True Strain'], 
                          **{header_90: avg_90_[header_90] for header_90 in avg_90_} }  
    
    elif (deg_0_file_red is None) and (deg_90_file_red is not None):
                    
        new_0_180_avg = {'Pattern-id':deg_0_file['Pattern-id'],
                         'Macro-True Strain':deg_0_file['Macro-True Strain'], 
                          **{header: avg_[header] for header in avg_} }


        new_90_270_avg = {'Pattern-id':deg_90_file_red['Pattern-id'],
                            'Macro-True Strain_90':deg_90_file_red['Macro-True Strain'], 
                         **{header_90: avg_90_[header_90] for header_90 in avg_90_} }
     
    elif (deg_0_file_red is not None) and (deg_90_file_red is None):
                     
         new_0_180_avg = {'Pattern-id':deg_0_file['Pattern-id'],
                         'Macro-True Strain':deg_0_file['Macro-True Strain'], 
                           **{header: avg_[header] for header in avg_} }

         new_90_270_avg = {'Pattern-id':deg_90_file['Pattern-id'],
                         'Macro-True Strain_90':deg_90_file['Macro-True Strain'], 
                          **{header_90: avg_90_[header_90] for header_90 in avg_90_} }
    
    else:
        
        new_0_180_avg = {'Pattern-id':deg_0_file['Pattern-id'],
                            'Macro-True Strain':deg_0_file['Macro-True Strain'], 
                          **{header: avg_[header] for header in avg_} }
        
        new_90_270_avg = {'Pattern-id':deg_90_file['Pattern-id'],
                            'Macro-True Strain_90':deg_90_file['Macro-True Strain'], 
                          **{header_90: avg_90_[header_90] for header_90 in avg_90_} }
    
    with pd.ExcelWriter(file_path + 'Average_microstrain_data.xlsx', engine='openpyxl') as out_file:
        pd.DataFrame(new_0_180_avg).to_excel(out_file, sheet_name='avg_0_180')
        pd.DataFrame(new_90_270_avg).to_excel(out_file, sheet_name='avg_90_270')
    
    print('Average file written in the path:', file_path)
        
    return 
        
oper_direc_file(file_path, sheet_name, alpha_beta_planes_dic)