#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a script to determine phase-level and material global level, 
single elastic constants (SEC) from microstrain and macroscopic stress vs strain 
data obtained from an in-situ high energy Xray diffraction (HEXRD).

This algorithm was developed by Muritala Oluwaseun Arowolo during his
PhD work at the Université de Toulouse, Toulouse France within two CNRS laboratories
(a) CIRIMAT - Toulouse INP
(b) ICA - Ecole des Mines Albi

The code handles SEC of HCP and BCC crystals. It can be modified to handle other
crystals.
"""
import os
import time
import numpy as np
import matplotlib.pyplot as plt
# import pygad
from scipy import stats
import pandas as pd
from numpy.linalg import inv
from numpy import linalg as LA
import random
from deap import creator, base, tools, algorithms
import sys
import csv


random.seed(680)
start = time.time()


file_path = "E:/New_PhD_Dell_19-03-2025/Result/New Petra/Ti6242-4mm-840C-103S-V2/2D-Integration/0-Degree/"


stress_strain_file = file_path + "SEC_microstrain_macrostress_macrostrain_data.xlsx"

micro_strain_sheet = 'Microstrain_FIT_Elastic'  
stress_strain_sheet = 'Stress-exp-elastic' 


'''Provide limit of experimental Stress vs strain to be 
    considered in the elastic region 
'''
lower_ind = 0
upper_ind = 80


'''Provide the phase SECs upper and lower bounds.
   ALPHA_BOUNDS are bounds of HCP phase SEC: C11,C33,C12,C13,C44 GPa
'''
ALPHA_BOUNDS = [(58.8, 168), (66.85, 191), (33.25, 95),
                 (24, 69), (15.8, 48)]  # C11,C33,C12,C13,C44 GPa ####


#BETA_BOUNDS are bounds of BCC phase SEC: C11,C12,C44  GPa
BETA_BOUNDS = [(26.9, 84), (28.5, 110), (15.25, 55)]   # C11,C12,C44  GPa 

#HCP c over a
c_over_a = 4.7374/2.9622 # 

#Volume fraction
vf_A = 0.82 ###Alpha - HCP
vf_B = 0.18 ###Beta - BCC

'''The following planes are the ones considered in the for 
  HCP and BCC crystals.
  It can be adjusted to handle more planes.
'''
a_1 = ["Alpha -{10-10}", "Alpha -{10-10}_90"]
a_2 = ["Alpha -{0002}", "Alpha -{0002}_90"]
a_3 = ["Alpha -{10-11}", "Alpha – {10-11}_90"]
a_4 = ["Alpha -{10-12}", "Alpha – {10-12}_90"]
a_5 = ["Alpha -{11-20}", "Alpha -{11-20}_90"]

b_1 = ["Beta -{110}", "Beta -{110}_90"]
b_2 = ["Beta -{200}", "Beta -{200}_90"]


alpha_planes_dic = {'a_1': a_1, 'a_2': a_2, 'a_3': a_3,
                    'a_4': a_4, 'a_5': a_5}
beta_planes_dic = {'b_1': b_1, 'b_2': b_2}


# microstrain
raw_micro_file = pd.read_excel(stress_strain_file, sheet_name=micro_strain_sheet.strip())

raw_micro = raw_micro_file.replace(
    ['', '-', '--', 'NA', 'N/A', ' ', 'nan', '#N/A', '#VALUE!', 'None'], np.nan)
micro_strain_df = raw_micro.astype(float)

#Data clean-up 
raw_stress = pd.read_excel(stress_strain_file, sheet_name=stress_strain_sheet.strip())

raw_stress = raw_stress.replace(
    ['', '-', '--', 'NA', 'N/A', ' ', 'nan', '#N/A',
        '#VALUE!', 'None', b'', b'-', b'--'],
    np.nan)

glb_stress = raw_stress.astype(float)

glb_stress_df = pd.DataFrame(glb_stress)


col_stress_exp_id = 'True stress' 
col_Glb_exp_strain_id = 'Strain_comp' 

glb_stress_exp = glb_stress_df[col_stress_exp_id]
glb_strain_exp = glb_stress_df[col_Glb_exp_strain_id]


# MICROSTRAIN Dataframe
micro_strain_df = pd.DataFrame(micro_strain_df) 

micro_strain_data = micro_strain_df



def fit_stress(strain_exp, stress_exp, micro_strain_data):
 
    linear_model = stats.linregress(strain_exp[lower_ind:upper_ind], 
                    stress_exp[lower_ind:upper_ind]-stress_exp[0])
    
    strain_pattern = micro_strain_data['Macro-True Strain']

    stress_at_strain_eqv_micro = linear_model.slope*(strain_pattern) + linear_model.intercept
    
    slp = (linear_model.slope/1000)

    plt.plot(strain_pattern, stress_at_strain_eqv_micro, label=f'Smooth - slope({slp:.2f}GPa)', color='black')
    plt.scatter(strain_exp[lower_ind:upper_ind], stress_exp[lower_ind:upper_ind], label='exp', alpha=0.2)
    plt.legend()
    plt.ylabel('Stress')
    plt.xlabel('Strain')
    plt.show

    return stress_at_strain_eqv_micro 



def stress_tensor_3by3_to_voigt(tensor):

    voigt_form = np.zeros(6)

    voigt_form[0] = tensor[0, 0]
    voigt_form[1] = tensor[1, 1]
    voigt_form[2] = tensor[2, 2]
    voigt_form[3] = 0.5 * (tensor[1, 2] + tensor[2, 1])
    voigt_form[4] = 0.5 * (tensor[0, 2] + tensor[2, 0])
    voigt_form[5] = 0.5*(tensor[0, 1] + tensor[1, 0])

    return voigt_form

def strain_tensor_3by3_to_voigt(tensor):

    strain_voigt_form = np.zeros(6)

    strain_voigt_form[0] = tensor[0, 0]
    strain_voigt_form[1] = tensor[1, 1]
    strain_voigt_form[2] = tensor[2, 2]
    strain_voigt_form[3] = 2*tensor[1, 2]
    strain_voigt_form[4] = 2*tensor[0, 2]
    strain_voigt_form[5] = 2*tensor[0, 1]

    return strain_voigt_form

def stiff_tensor(C11:float, C12, C44, C33=None ,C13=None):
    
    tensor = np.zeros((6,6))
    if C33 is None:
        tensor[0, 0] = tensor[1, 1] = tensor[2, 2] = C11
        tensor[0, 1] = tensor[1, 0] = C12
        tensor[0, 2] = tensor[2, 0] = C12
        tensor[1, 2] = tensor[2, 1] = C12
        tensor[3, 3] = tensor[4, 4] = tensor[5, 5] = C44

    else:
        tensor[0, 0] = tensor[1, 1] = C11
        tensor[2, 2] = C33 #stiffness[0]
        tensor[0, 1] = tensor[1, 0] = C12
        tensor[0, 2] = tensor[2, 0] = C13
        tensor[1, 2] = tensor[2, 1] = C13
        tensor[3, 3] = tensor[4, 4] = C44
        tensor[5, 5] = 0.5 * (C11 - C12)
      
    return tensor


def return_val(tensor):
    
    C11 = tensor[0, 0]
    C12 = tensor[0, 1]
    C44 = tensor[3, 3] 
    C33 = tensor[2, 2]
    C13 = tensor[0, 2] 

    return C11, C12, C44, C33,C13



def stress_exp_comp(stress):

    stress_tensor = np.asmatrix(np.zeros((3, 3)))
    stress = float(stress)

    stress_tensor[0, 0] = stress

    return stress_tensor


def von_mises(inpt):

    a = 0.5*((inpt[0,0] - inpt[1,1])**2 + (inpt[1,1] - inpt[2,2])**2  + (inpt[2,2] - inpt[0,0])**2)

    b = 3* ((0.5*(inpt[0,1]+inpt[1,0]))**2 + (0.5*(inpt[1,2]+inpt[2,1]))**2 +
            (0.5*(inpt[0,2]+inpt[2,0]))**2)

    vm_eq = np.sqrt(a + b)

    return vm_eq

def trial_hkl_sec(inp):
    h, k, l = inp[0], inp[1], inp[2]

    if (h, k, l) == (0, 0, 0):
        raise ValueError("Invalid plane (0,0,0)")
        sys.exit()
    check = np.nonzero(inp)[0]
    
    if np.all(inp != 0):
        u = k
        v = h
        w = -(h*u + k*v)/l

    elif (list(check) == [0,1]):
        u = -k
        v = h
        w = 0
    elif (list(check) == [1,2]):
        u = 0
        v = -l
        w = k
    elif (list(check) == [0,2]):
        u = -l
        v = 0
        w = h       
        
    elif len(check)==1:
        if (list(check) == [0]):
            u = 0
            v = -h
            w = 0
        elif (list(check) == [1]):
            u = -k
            v = 0
            w = 0
        else:
            u = -l
            v = 0
            w = 0
    
    else:  
        print('Unknown hkl, so uvw that satisfies Weiss law cannot be found')
        sys.exit(0)

    return u, v, w   


def strain_rotation_tensor(plane, ca):
    root_2 = np.sqrt(2)
    elm = 1.0 / root_2
    
    tr_mat = np.array([
        [1.0,  0.0, 0.0],
        [-1/2,  np.sqrt(3)/2,  0.0],
        [0.0, 0.0,  ca]
        ])

    if plane.strip() == 'a_1':        # Alpha -{10-10} prismatic
    
        eqv_plane = [1, 0, 0]         # Alpha -{100}
        
        Tr_plane = np.matmul(tr_mat , np.transpose(eqv_plane))
        mag_plane = LA.norm(Tr_plane)
  
        dir_Tr_plane = Tr_plane/mag_plane
       
        h,k,l = dir_Tr_plane
        u,v,w = trial_hkl_sec(dir_Tr_plane)
        p,q, r = np.cross([h,k,l]  ,[u,v,w])
        
        R = np.array([
            [h, u, p],
            [k, v, q],
            [l, w, r]])

    elif plane == "a_2":      # Alpha -{0002} basal
        
        eqv_plane = [0, 0, 2]         # Alpha -{002}
        
        Tr_plane = np.matmul(tr_mat , np.transpose(eqv_plane))

        mag_plane = LA.norm(Tr_plane)

        dir_Tr_plane = Tr_plane/mag_plane
       
        h,k,l = dir_Tr_plane

        u,v,w = trial_hkl_sec(dir_Tr_plane)

        p,q, r = np.cross([h,k,l]  ,[u,v,w])
        
        R = np.array([
            [h, u, p],
            [k, v, q],
            [l, w, r]])

    elif plane == 'a_3':      # Alpha -{10-11} pyramidal 
    
        eqv_plane = [1, 0, 1]         # Alpha -{10-11}
        
        Tr_plane = np.matmul(tr_mat , np.transpose(eqv_plane))
        mag_plane = LA.norm(Tr_plane)
        
        dir_Tr_plane = Tr_plane/mag_plane
       
        h,k,l = dir_Tr_plane
        u,v,w = trial_hkl_sec(dir_Tr_plane)
        p,q, r = np.cross([h,k,l]  ,[u,v,w])
        
        R = np.array([
            [h, u, p],
            [k, v, q],
            [l, w, r]])
        #print('out', R)
    elif plane == 'a_4':      # Alpha -{10-12} pyramidal 
    
        eqv_plane = [1, 0, 2]         # Alpha -{10-12}
        
        Tr_plane = np.matmul(tr_mat , np.transpose(eqv_plane))
        mag_plane = LA.norm(Tr_plane)

        dir_Tr_plane = Tr_plane/mag_plane
       
        h,k,l = dir_Tr_plane
        u,v,w = trial_hkl_sec(dir_Tr_plane)
        p,q, r = np.cross([h,k,l]  ,[u,v,w])
        
        R = np.array([
            [h, u, p],
            [k, v, q],
            [l, w, r]])

    elif plane == 'a_5':      # Alpha -{11-20} prismatic 
        
        eqv_plane = [1, 1, 0]         # Alpha -{11-20}
        
        Tr_plane = np.matmul(tr_mat , np.transpose(eqv_plane))
        mag_plane = LA.norm(Tr_plane)
        #print(mag_plane)
        dir_Tr_plane = Tr_plane/mag_plane
       
        h,k,l = dir_Tr_plane
        u,v,w = trial_hkl_sec(dir_Tr_plane)
        p,q, r = np.cross([h,k,l]  ,[u,v,w])
        
        R = np.array([
            [h, u, p],
            [k, v, q],
            [l, w, r]])

    elif plane == 'b_1':      # Beta -{110}
        R = np.array([
            [elm, elm, 0.0],
            [elm, -elm, 0.0],
            [0.0, 0.0, 1.0]])
    
    elif plane == 'b_2':      # Beta -{200}
        R = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]])

    else:
        print(f'The plane is unknown or there is no matching rotation matrix for {plane} plane provided')
        R = np.eye(3)

    return R


def strain_tensor(s_11, s_22):

    strain_tensor = np.zeros((3, 3))

    strain_tensor[0, 0] = s_11
    strain_tensor[1, 1] = strain_tensor[2, 2] = s_22

    return strain_tensor


def get_plane_data(strain_data, plane):
    if plane is None:
        return strain_data.iloc[:, 0]
    try:
        return strain_data[plane]
    except:
        print(f"Plane '{plane}' not found using first column")
        return strain_data.iloc[:, 0]


def set_up_stiffness_matrix(phase, stiffness_list):
    C = np.zeros((6, 6))
    if phase.lower() == 'hcp':
        C11, C33, C12, C13, C44 = stiffness_list
        C66 = (C11 - C12)/2.0

        if C66 < 1e-6:
            C66 = 1e-6

        C[0, 0] = C[1, 1] = C11
        C[2, 2] = C33
        C[0, 1] = C[1, 0] = C12
        C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = C13
        C[3, 3] = C[4, 4] = C44
        C[5, 5] = C66

    elif phase.lower() == 'bcc':
        C11, C12, C44 = stiffness_list
        C[0:3, 0:3] = C12
        np.fill_diagonal(C[0:3, 0:3], C11)
        C[3, 3] = C[4, 4] = C[5, 5] = C44

    return C * 1000.0


def strain_transformation(strain_hkl, rotat_matrix):

    inv_rot_mat = np.linalg.inv(rotat_matrix)
    transf_strain = inv_rot_mat @ strain_hkl @ inv_rot_mat.T

    return transf_strain


def voigt_to_tensor(voigt_vector):

    v = np.asarray(voigt_vector).ravel()

    tensor = np.zeros((3, 3))
    tensor[0, 0] = v[0]   
    tensor[1, 1] = v[1]   
    tensor[2, 2] = v[2]   
    tensor[1, 2] = tensor[2, 1] = v[3]  
    tensor[0, 2] = tensor[2, 0] = v[4]  
    tensor[0, 1] = tensor[1, 0] = v[5]  

    return tensor


def target_stress(stress_exp_data, col_stress_exp_id, col_Glb_exp_strain_id,
                  micro_strain_data):

    stress_glb = get_plane_data(stress_exp_data, col_stress_exp_id)

    strain_glb = get_plane_data(stress_exp_data, col_Glb_exp_strain_id)

    smooth_stress = fit_stress(strain_glb, stress_glb, micro_strain_data)

    return np.array(smooth_stress), strain_glb



def voigt_model_stress(micro_strain_data, vf_A, vf_B, c_over_a,
                       alpha_stiffness_list, beta_stiffness_list,
                       alpha_planes_dic, beta_planes_dic):

    predicted_vght_von_mises = []

    C_alpha = set_up_stiffness_matrix('hcp', alpha_stiffness_list)
    C_beta = set_up_stiffness_matrix('bcc', beta_stiffness_list)

    n_steps = len(micro_strain_data)

    for i in range(n_steps):
        cum_alpha_strain = np.zeros((3, 3))
        cum_beta_strain = np.zeros((3, 3))

        for plane_key in alpha_planes_dic.keys():
            plane_0d = alpha_planes_dic[plane_key][0]
            plane_90d = alpha_planes_dic[plane_key][1]
            R_alpha = strain_rotation_tensor(plane_key, c_over_a)

            eps_11 = get_plane_data(micro_strain_data, plane_0d).iloc[i] ###Loading directon
            eps_22 = get_plane_data(micro_strain_data, plane_90d).iloc[i] ###Transverse directon


            eps_hkl = strain_tensor(eps_11, eps_22)

            eps_lab_alpha = strain_transformation(eps_hkl, R_alpha)
            cum_alpha_strain += eps_lab_alpha

        for plane_key in beta_planes_dic.keys():

            plane_0d = beta_planes_dic[plane_key][0]
            plane_90d = beta_planes_dic[plane_key][1]
            R_beta = strain_rotation_tensor(plane_key, c_over_a)

            eps_11 = get_plane_data(micro_strain_data, plane_0d).iloc[i]   ###Loading directon
            eps_22 = get_plane_data(micro_strain_data, plane_90d).iloc[i]  ###Transverse directon

            eps_hkl = strain_tensor(eps_11, eps_22)
            eps_lab_beta = strain_transformation(eps_hkl, R_beta)
            cum_beta_strain += eps_lab_beta

        eps_alpha_voigt = strain_tensor_3by3_to_voigt(cum_alpha_strain)
        eps_beta_voigt = strain_tensor_3by3_to_voigt(cum_beta_strain)

        sigma_alpha_voigt = C_alpha @ eps_alpha_voigt
        sigma_beta_voigt = C_beta  @ eps_beta_voigt

        sigma_total_voigt = vf_A * sigma_alpha_voigt + vf_B * sigma_beta_voigt
        sigma_tensor = voigt_to_tensor(sigma_total_voigt)

        predicted_vght_von_mises.append(von_mises(sigma_tensor))

    return predicted_vght_von_mises, np.array(micro_strain_data['Macro-True Strain'])


def is_stable_hcp(C11, C33, C12, C13, C44):
    """
    Necessary and sufficient elastic stability conditions for HCP
    Reference: Mouhat & Coudert, Phys. Rev. B 90, 224104 (2014)
    """

    C66 = 0.5 * (C11 - C12)
    cond1 = C11 > abs(C12)                    
    cond2 = C44 > 0                            
    cond3 = C66 > 0
    cond4 = 2 * C13**2 < C33 * (C11 + C12)
    cond5 = C33 > C11

    return cond1 and cond2 and cond3 and cond4 and cond5

def is_stable_bcc(C11, C12, C44):
    """
    Necessary and sufficient elastic stability conditions for HCP
    Reference: Mouhat & Coudert, Phys. Rev. B 90, 224104 (2014)
    """
    cond1 = (C11 - C12) > 0      
    cond2 = (C11 + 2*C12) > 0    
    cond3 = C44 > 0
    cond4 = C11 > C44 

    return cond1 and cond2 and cond3 and cond4

def limit_bound(ALPHA_BOUNDS, BETA_BOUNDS, individual):
    
    C11_A, C33_A, C12_A, C13_A, C44_A, C11_B, C12_B, C44_B = individual
    
    #ALPHA UPPER LIMIT
    tot_cond = [
    ALPHA_BOUNDS[0][0] <= C11_A < ALPHA_BOUNDS[0][1], 
    ALPHA_BOUNDS[1][0] <= C33_A < ALPHA_BOUNDS[1][1],
    ALPHA_BOUNDS[2][0] <= C12_A < ALPHA_BOUNDS[2][1],
    ALPHA_BOUNDS[3][0] <= C13_A < ALPHA_BOUNDS[3][1], 
    ALPHA_BOUNDS[4][0] <= C44_A < ALPHA_BOUNDS[4][1],

    #BETA UPPER LIMIT
    BETA_BOUNDS[0][0]<= C11_B < BETA_BOUNDS[0][1],
    BETA_BOUNDS[1][0]<= C12_B < BETA_BOUNDS[1][1],
    BETA_BOUNDS[2][0]<= C44_B < BETA_BOUNDS[2][1]

    ]

    return all(tot_cond)

def global_elastic_B_v_E(stiffness):

    C11, C33, C12, C13, C44 = stiffness[:5]
    hcp_stiffness = stiff_tensor(C11, C12, C44, C33, C13)

    C11_, C12_, C44_ = stiffness[5:]
    bcc_stiffness = stiff_tensor(C11_, C12_, C44_)

    C_Glob = vf_A * hcp_stiffness + vf_B * bcc_stiffness

    C11_G, C12_G, C44_G, C33_G, C13_G = return_val(C_Glob)

    E = ((C33_G*(C11_G + C12_G) - 2*C13_G**2) * (C11_G - C12_G)) / \
        ((C33_G*C11_G) - C13_G**2)      ###Guler et al. 2021
    
    K = ((2*C11_G + C33_G) + 2*(C12_G + 2*C13_G))/9    ### Hill 1952
    
    C66_G = (C11_G - C12_G)/2
    
    G = ((2*C11_G + C33_G) - (C12_G + 2*C13_G) +  3*(2*C44_G + C66_G))/15 ### Hill 1952
    
    v = 0.5*(1 - (3*G)/(3*K +G))        ### Hill 1952
    
    hcp_dic  = {'HCP_C11': C11, 'HCP_C12': C12, 'HCP_C44': C44, 'HCP_C33': C33, 'HCP_C13': C13,}
    
    bcc_dic  = {'BCC_C11': C11_, 'BCC_C12': C12_, 'BCC_C44': C44_}
    
    res_Glob = {
        "C11_G": round(C11_G, 2),
        "C12_G": round(C12_G, 2),
        "C44_G": round(C44_G, 2),
        "C33_G": round(C33_G, 2),
        "C13_G": round(C13_G, 2),
        "E": round(E, 2),
        "G": round(G, 2),
        "K": round(K, 2),
        "v": round(v, 2)     
    }
    
    print(res_Glob)
    
    sav_path = os.path.join(file_path, 'Single_elastic_constant_output.txt')
    with open(sav_path, 'w', newline='') as file:
        wrt = csv.writer(file, delimiter='\t')
    
        wrt.writerow(['--- HCP Elastic Constants ---'])
        wrt.writerows(hcp_dic.items())
    
        wrt.writerow([])
        wrt.writerow(['--- BCC Elastic Constants ---'])
        wrt.writerows(bcc_dic.items())
    
        wrt.writerow([])
        wrt.writerow(['--- Global Elastic Constants ---'])
        wrt.writerows(res_Glob.items())
    
    return 


# ==============================================================
# GENETIC ALGORITHM WITH DEAP 
# ==============================================================

target_vm, strain_glb = target_stress(glb_stress_df, col_stress_exp_id, col_Glb_exp_strain_id,
                                      micro_strain_df)

####Initialization
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Register attributes
attrs = ["C11_A", "C33_A", "C12_A", "C13_A",
         "C44_A", "C11_B", "C12_B", "C44_B"]
all_bounds = ALPHA_BOUNDS + BETA_BOUNDS


for name, (lo, hi) in zip(attrs, all_bounds):
    toolbox.register(name, random.uniform, lo, hi)

toolbox.register("individual", tools.initCycle, creator.Individual,
                 (toolbox.C11_A, toolbox.C33_A, toolbox.C12_A, toolbox.C13_A,
                  toolbox.C44_A, toolbox.C11_B, toolbox.C12_B, toolbox.C44_B), n=1)

toolbox.register("population", tools.initRepeat, list, toolbox.individual)


def evaluate_stiffness(individual):
    C11_A, C33_A, C12_A, C13_A, C44_A, C11_B, C12_B, C44_B = individual
    
    if not limit_bound(ALPHA_BOUNDS, BETA_BOUNDS, individual):
        return (1e15,)
    if not is_stable_hcp(C11_A, C33_A, C12_A, C13_A, C44_A):
        return (1e15,)
    if not is_stable_bcc(C11_B, C12_B, C44_B):
        return (1e15,)


    alpha_stiffness = [C11_A, C33_A, C12_A, C13_A, C44_A]
    beta_stiffness = [C11_B, C12_B, C44_B]

    try:
        predicted, _  = voigt_model_stress(
            micro_strain_data, vf_A, vf_B, c_over_a,
            alpha_stiffness, beta_stiffness,
            alpha_planes_dic, beta_planes_dic)

        mse = np.mean((predicted - target_vm)**2)

        return (mse,) 

    except Exception as e:
        print(f"CRASH - {e}")
        return (1e15,)


toolbox.register("evaluate", evaluate_stiffness)
toolbox.register("select", tools.selTournament, tournsize=5)
toolbox.register("mate", tools.cxBlend, alpha=0.5)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=12, indpb=0.25)

pop = toolbox.population(n= 450)
hof = tools.HallOfFame(1)
stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
stats.register("min", np.min)
stats.register("avg", np.mean)

pop, log = algorithms.eaSimple(
    pop, toolbox,
    cxpb=0.7, mutpb=0.7,
    ngen=350,
    stats=stats,
    halloffame=hof,
    verbose=False)

# ==============================================================
# 6. BEST RESULT
# ==============================================================

best = hof[0]
print("\n" + "="*70)
print("OPTIMIZATION COMPLETE — BEST ELASTIC CONSTANTS FOUND")
print("="*70)
print(f"α-phase (HCP):")
print(f"  C11 = {best[0]:6.1f} GPa,  C33 = {best[1]:6.1f} GPa")
print(f"  C12 = {best[2]:6.1f} GPa,  C13 = {best[3]:6.1f} GPa,  C44 = {best[4]:6.1f} GPa")
print(f"β-phase (BCC):")
print(f"  C11 = {best[5]:6.1f} GPa,  C12 = {best[6]:6.1f} GPa,  C44 = {best[7]:6.1f} GPa")
print(f"\nFinal MSE = {best.fitness.values[0]:.3f} MPa²")
print("="*70)
print('target_vm', target_vm)


global_elastic_B_v_E(best)

alpha_stiffness = list(best[:5]) #[C11_A, C33_A, C12_A, C13_A, C44_A]
beta_stiffness = list(best[5:]) #[C11_B, C12_B, C44_B]

predicted, strain__  = voigt_model_stress(
    micro_strain_data, vf_A, vf_B, c_over_a,
    alpha_stiffness, beta_stiffness,
    alpha_planes_dic, beta_planes_dic)
print(predicted)

plt.figure()
plt.plot(strain__, predicted, label='Sim')
plt.ylabel('Stress')
plt.xlabel('Strain')
plt.scatter(strain__, target_vm, label='Target')
plt.legend()
plt.show()
end = time.time()
print(f'Total computation time, {(end - start):.2f} s')

