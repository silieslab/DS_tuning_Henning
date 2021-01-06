#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

@author: burakgur

Required functions to analyze T4 T5 delays


"""
import numpy as np
import matplotlib.pyplot as plt

def run_matplotlib_params():
    plt.style.use('default')
    plt.style.use('seaborn-talk')
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)    
    plt.rcParams["axes.titlesize"] = 'medium'
    plt.rcParams["axes.labelsize"] = 'small'
    plt.rcParams["axes.labelweight"] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams["legend.fontsize"] = 'small'
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["figure.titleweight"] = 'bold'
    plt.rcParams["figure.titlesize"] = 'medium'
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.fontsize'] = 'x-small'
    plt.rcParams['legend.loc'] = 'upper right'
    
    c_dict = {}
    c_dict['dark_gray'] = np.array([77,77,77]).astype(float)/255
    c_dict['light_gray'] = np.array([186,186,186]).astype(float)/255
    c_dict['green1'] = np.array([102,166,30]).astype(float)/255
    c_dict['green2']=np.array([179,226,205]).astype(float)/255
    c_dict['green3'] = np.array([27,158,119]).astype(float)/255
    c_dict['orange']  = np.array([201,102,47]).astype(float)/255
    c_dict['red']  = np.array([228,26,28]).astype(float)/255
    c_dict['magenta']  = np.array([231,41,138]).astype(float)/255
    c_dict['purple']  = np.array([117,112,179]).astype(float)/255
    c_dict['yellow'] = np.array([255,255,51]).astype(float)/255
    c_dict['brown'] = np.array([166,86,40]).astype(float)/255

    
    c = []
    c.append(c_dict['dark_gray'])
    c.append(c_dict['light_gray'])
    c.append(c_dict['green1']) # Green
    c.append(c_dict['orange']) # Orange
    c.append(c_dict['red']) # Red
    c.append(c_dict['magenta']) # magenta
    c.append(c_dict['purple'])# purple
    c.append(c_dict['green2']) # Green
    c.append(c_dict['yellow']) # Yellow
    c.append(c_dict['brown']) # Brown
    
    
    return c, c_dict

def apply_threshold_df(threshold_dict, df):
    
    if threshold_dict is None:
        print('No threshold used.')
        return df
    
    pass_bool = np.ones((1,len(df)))
    
    for key, value in threshold_dict.items():

        pass_bool = pass_bool * np.array((df[key] > value))
        
    threshold_df = df[pass_bool.astype(bool)[0]]
   
    return threshold_df

def threshold_ROIs(rois, threshold_dict):
    """ Thresholds given ROIs and returns the ones passing the threshold.

    Parameters
    ==========
    rois : list
        A list of ROI_bg instances.
        
    threshold_dict: dict
        A dictionary with desired ROI_bg property names that will be 
        thresholded as keys and the corresponding threshold values as values. 
    
    Returns
    =======
    
    thresholded_rois : list 
        A list containing instances of ROI_bg which pass the thresholding step.
    """
    # If there is no threshold
    if threshold_dict is None:
        print('No threshold used.')
        return rois
    vars_to_threshold = threshold_dict.keys()
    
    roi_data_dict = data_to_list(rois, vars_to_threshold)
    
    pass_bool = np.ones((1,len(rois)))
    
    for key, value in threshold_dict.items():
        
        if type(value) == tuple:
            if value[0] == 'b':
                pass_bool = \
                    pass_bool * (np.array(roi_data_dict[key]).flatten() > value[1])
                
            elif value[0] == 's':
                pass_bool = \
                    pass_bool * (np.array(roi_data_dict[key]).flatten() < value[1])
            else:
                raise TypeError("Tuple first value not understood: should be 'b' for bigger than or 's' for smaller than")
                
        else:
            pass_bool = pass_bool * (np.array(roi_data_dict[key]).flatten() > value)
    
    pass_indices = np.where(pass_bool)[1]
    
    thresholded_rois = []
    for idx in pass_indices:
        thresholded_rois.append(rois[idx])
    
    return thresholded_rois

def data_to_list(rois, data_name_list):
    """ Generates a dictionary with desired variables from ROIs.

    Parameters
    ==========
    rois : list
        A list of ROI_bg instances.
        
    data_name_list: list
        A list of strings with desired variable names. The variables should be 
        written as defined in the ROI_bg class. 
        
    Returns
    =======
    
    roi_data_dict : dictionary 
        A dictionary with keys as desired data variable names and values as
        list of data.
    """   
    class my_dictionary(dict):  
  
        # __init__ function  
        def __init__(self):  
            self = dict()  
              
        # Function to add key:value  
        def add(self, key, value):  
            self[key] = value  
    
    roi_data_dict = my_dictionary()
    
    # Generate an empty dictionary
    for key in data_name_list:
        roi_data_dict.add(key, [])
    
    # Loop through ROIs and get the desired data            
    for iROI, roi in enumerate(rois):
        for key, value in roi_data_dict.items(): 
            if key in roi.__dict__.keys():
                value.append(roi.__dict__[key])
            else:
                value.append(np.nan)
    return roi_data_dict

def generate_time_delay_profile_2Dedges(rois):
    
    for roi in rois:
        epochDur= roi.stim_info['epochs_duration']
        max_epoch = roi.max_resp_idx
        roi.edge_start_loc = roi.stim_info['input_data']['Stimulus.stimtrans.amp'][max_epoch]
        roi.edge_speed = roi.stim_info['input_data']['Stimulus.stimtrans.mean'][max_epoch]
        half_dur_frames = int((round(roi.imaging_info['frame_rate'] * epochDur[max_epoch]))/2)
        trace = roi.resp_trace_all_epochs[roi.max_resp_idx]
        OFF_resp = trace[:half_dur_frames]
        ON_resp = trace[half_dur_frames:]
        if roi.CS =='ON':
            roi.two_d_edge_profile = ON_resp
        elif roi.CS =='OFF':
            roi.two_d_edge_profile = OFF_resp
    return rois

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def fit_1d_gauss(data_x, data_y):
    from scipy.optimize import curve_fit
    p0 = [np.max(data_y), np.argmax(data_y), 1]
    coeff, pcov = curve_fit(gauss, data_x, data_y, p0=p0)
    fit_trace = gauss(data_x, *coeff)
    
    residuals = data_y- fit_trace
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((data_y-np.mean(data_y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return fit_trace, r_squared, coeff

def generate_time_delay_profile_combined(rois,screen_deg = 60):
    
    screen_coords = np.linspace(0, screen_deg, num=screen_deg, endpoint=True) # degree of the screen
    
    # Edge is presented in the full screen and not just 60 degrees of the visual field

    for roi in rois:
        
        roi_t_v_stripe = np.linspace(0, screen_deg, num=len(roi.max_resp_all_epochs[1:]), endpoint=True)
        
        
        diff = np.abs(int(roi.edge_start_loc)) - 40 -screen_deg/2
        start_frame = int(np.around((diff/float(roi.edge_speed)) * roi.imaging_info['frame_rate']))
        end_frame = start_frame + int(np.ceil((60/float(roi.edge_speed) * roi.imaging_info['frame_rate'])))
        roi_t_v_edge = \
            np.linspace(0, screen_deg, 
                        num=len(roi.two_d_edge_profile[start_frame:end_frame]),
                        endpoint=True)
        
        i_edge = np.interp(screen_coords, roi_t_v_edge, 
                           roi.two_d_edge_profile[start_frame:end_frame])
        roi.i_stripe_resp = np.interp(screen_coords, roi_t_v_stripe, 
                                  np.transpose(roi.max_resp_all_epochs[1:])[0])
        if roi.PD == 90: # Rotate the response if PD is 90   
            roi.i_edge_resp = i_edge[::-1]
        else:
            roi.i_edge_resp = i_edge
        try:
            fit_trace, r_squared, coeff = fit_1d_gauss(screen_coords, roi.i_edge_resp)
        except RuntimeError:
            print('Fit parameters not found... discarding {s}'.format(s=roi))
            roi.discard = True
            roi.edge_gauss_profile = None
            roi.edge_r_squared = None
            roi.edge_gauss_coeff = None
            roi.resp_delay_deg = None
            roi.resp_delay_fits_Rsq = None
            continue
        roi.edge_gauss_profile = fit_trace
        roi.edge_r_squared = r_squared
        roi.edge_gauss_coeff = coeff
        
        try:
            fit_trace, r_squared, coeff = fit_1d_gauss(screen_coords, roi.i_stripe_resp)
        except RuntimeError:
            print('Fit parameters not found... discarding {s}'.format(s=roi))
            roi.discard = True
            roi.stripe_gauss_profile = None
            roi.stripe_r_squared = None
            roi.stripe_gauss_coeff = None
            roi.resp_delay_deg = None
            roi.resp_delay_fits_Rsq = None
            continue
        roi.stripe_gauss_profile = fit_trace
        roi.stripe_r_squared = r_squared
        roi.stripe_gauss_coeff = coeff
        
        roi.resp_delay_deg = np.abs(np.argmax(roi.stripe_gauss_profile) - np.argmax(roi.edge_gauss_profile))
        roi.resp_delay_fits_Rsq = np.array([roi.edge_r_squared,roi.stripe_r_squared])
        roi.discard = False
        
    return rois