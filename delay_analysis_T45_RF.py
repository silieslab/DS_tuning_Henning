#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:08:09 2020

@author: burakgur

Analysis

"""
#%% Importing required packages
import cPickle
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import warnings

os.chdir('/Users/burakgur/Documents/GitHub/python_lab/2p_calcium_imaging')
import henning_helper as hh

# %% Setting the directories
initialDirectory = 'DS_tuning_Henning'
datasets_dir = os.path.join(initialDirectory, 'Data','delay_calculation_data')

# %% Load datasets and desired variables
datasets_to_load = ['191209bg_fly1-TSeries-002_transfer_sima_STICA.pickle',
                    '191209bg_fly1-TSeries-003_transfer_sima_STICA.pickle',
                    '191209bg_fly2-TSeries-004_transfer_sima_STICA.pickle',
                    '191209bg_fly2-TSeries-005_transfer_sima_STICA.pickle',
                    '191210bg_fly1-TSeries-12102019-0944-009_transfer_sima_STICA.pickle',
                    '191210bg_fly1-TSeries-12102019-0944-010_transfer_sima_STICA.pickle',
                    '191210bg_fly2-TSeries-12102019-0944-002_transfer_sima_STICA.pickle',
                    '191210bg_fly2-TSeries-12102019-0944-003_transfer_sima_STICA.pickle',
                    '191210bg_fly2-TSeries-12102019-0944-009_transfer_sima_STICA.pickle',
                    '191210bg_fly2-TSeries-12102019-0944-010_transfer_sima_STICA.pickle',
                    '191212bg_fly2-TSeries-002_transfer_sima_STICA.pickle',
                    '191212bg_fly2-TSeries-004_transfer_sima_STICA.pickle',
                    '191212bg_fly3-TSeries-004_transfer_sima_STICA.pickle',
                    '191212bg_fly3-TSeries-005_transfer_sima_STICA.pickle',
                    '191213bg_fly2-TSeries-12132019-0909-005_transfer_sima_STICA.pickle',
                    '191213bg_fly2-TSeries-12132019-0909-006_transfer_sima_STICA.pickle',
                    '191213bg_fly3-TSeries-12132019-0909-004_transfer_sima_STICA.pickle',
                    '191216bg_fly3-TSeries-004_transfer_sima_STICA.pickle',
                    '191216bg_fly3-TSeries-005_transfer_sima_STICA.pickle',
                    ]
properties = ['Delay', 'Rsq', 'CSI','Reliab','CS']
combined_df = pd.DataFrame(columns=properties)
all_rois = []
for idataset, dataset in enumerate(datasets_to_load):

    # Load the ROIs from .pickle files
    load_path = os.path.join(datasets_dir, dataset)
    load_path = open(load_path, 'rb')
    workspace = cPickle.load(load_path)
    curr_rois = workspace['final_rois']

    # Find the delay between moving and static stimulus 
    curr_rois=hh.generate_time_delay_profile_combined(curr_rois)


    # Threshold based on reliability of each ROI. 
    # Reliability is the correlation of responses for a given epoch in different trials. 
    # High correlation means that an ROI was responding reliably.
    curr_rois = hh.threshold_ROIs(curr_rois, {'reliability' : 0.55})
    if len(curr_rois) < 1:
        warnings.warn('{d} contains no ROIs, skipping'.format(d=dataset))
        continue

    data_to_extract = ['SNR', 'CSI','reliability','resp_delay_fits_Rsq','PD',
                       'resp_delay_deg','CS']
    roi_data = hh.data_to_list(curr_rois, data_to_extract)
    
    # There are two fits (moving and static epochs) so take the minimum R-squared value for a conservative estimate for each ROI.
    # An ROI can have None as delay or Rsq if it is not responding
    rsq = np.array(map(np.nanmin,roi_data['resp_delay_fits_Rsq'])) 
    rsq[rsq==None] = np.nan
    rsq[rsq<0] = np.nan
    
    delay = np.array(map(np.nanmin,roi_data['resp_delay_deg']))
    if len(delay[delay==None]):
        delay[delay==None] = np.nan
    delay = delay.astype(float)

    pref_dir = roi_data['PD']
    pd_S = map(str,list(map(int,pref_dir)))
    
    # Collect everything in a common pandas dataframe
    df_c = {}
    df_c['Delay'] = delay
    df_c['Rsq'] = rsq
    df_c['PD'] = pd_S
    df_c['SNR'] = roi_data['SNR']
    df_c['CSI'] = roi_data['CSI']
    df_c['Reliab'] = roi_data['reliability']
    df_c['CS'] = roi_data['CS']
    df_c['flyID'] = np.tile(curr_rois[0].experiment_info['FlyID'],len(curr_rois))
    
    df = pd.DataFrame.from_dict(df_c) 
    rois_df = pd.DataFrame.from_dict(df)
    combined_df = combined_df.append(rois_df, ignore_index=True, sort=False)
    
    all_rois.append(curr_rois)
    
    print('{ds} successfully loaded\n'.format(ds=dataset))
all_rois = np.concatenate(all_rois)
#%% Thresholding to get rid of bad fits
rsq_t = 0.7
threshold_dict = {'Rsq': rsq_t}
threshold_df = hh.apply_threshold_df(threshold_dict, combined_df)

#%% Generate summary figure
plt.close('all')
_, colorss = hh.run_matplotlib_params()
colors = [colorss['magenta'], colorss['green1']]

fig = plt.figure(figsize=(14, 3))

grid = plt.GridSpec(1,3 ,wspace=0.3, hspace=1)

ax1=plt.subplot(grid[0,0])
ax2=plt.subplot(grid[0,1])
ax3=plt.subplot(grid[0,2])

sns.scatterplot(x='Delay',y='Rsq',s=20, linewidth=0,alpha=.3,
             ax=ax1,data=combined_df,hue='Reliab',palette='inferno')

# sns.scatterplot(x='Delay',y='Rsq',s=20, linewidth=0,alpha=.5,hue='Reliab',
             # ax=ax1,data=threshold_df)

plt.setp(ax1.get_legend().get_texts(), fontsize='6')
ax1.plot([0, 60],[rsq_t,rsq_t],'--',color = 'k',alpha=.8)
ax1.set_xlim((0,60))
ax1.set_ylabel('R$^2$')
ax1.set_xlabel('Delay ($^\circ$)')
ax1.set_title('Thresholded: {p}/{a} ROIs'.format(p = len(threshold_df), a=len(combined_df)))
for idx, cs in enumerate(np.unique(threshold_df['CS'])):
    curr_df = threshold_df[threshold_df['CS'] == cs]
    lb = '{cs} {f} ({roi}) '.format(cs=cs,f=len(np.unique(curr_df['flyID'])),
                                   roi=len(curr_df))
    sns.distplot(curr_df['Delay'],ax=ax2,kde=True,
                 kde_kws={"alpha": .7,'cumulative': False},
                 color=(colors[idx][0],colors[idx][1],colors[idx][2],1.0),
                     hist=False,bins=10,label=lb)
    
ax2.set_xlabel('Delay ($^\circ$)')
ax2.set_title('Edge delay from RF center')

sns.violinplot(x='CS',y='Delay',width=0.5,ax=ax3,linewidth=1.5,data=threshold_df,
               palette=colors,inner="quartile")
# ax3.set_ylim((0, 25))
ax3.set_xlabel('')
ax3.set_ylabel('Delay ($^\circ$)')
ax3.legend('')
mon = threshold_df[threshold_df['CS']=='ON']['Delay'].mean()
moff = threshold_df[threshold_df['CS']=='OFF']['Delay'].mean()
son = threshold_df[threshold_df['CS']=='ON']['Delay'].std()
soff = threshold_df[threshold_df['CS']=='OFF']['Delay'].std()
ax3.set_title('OFF: {moff}$^\circ$$\pm${soff}$^\circ$   ON: {mon}$^\circ$$\pm${son}$^\circ$'.format(moff=np.around(moff,1),
                                                           mon=np.around(mon,1),
                                                           soff=np.around(soff,1),
                                                           son=np.around(son,1)))
fig.tight_layout()



# %%
