The codes were written for the following paper:

MY Yim, X Cai & X-J Wang (2019) Transforming the choice outcome to an action plan in monkey lateral prefrontal cortex: a neural circuit model.

The MATLAB codes are written in MATLAB R2014b by Man Yi Yim (manyi.yim@gmail.com) and Xinying Cai (xinying.cai@nyu.edu). The Python codes are written in Python 2.7 by Man Yi Yim. The codes were finalized on May 24, 2019.

### Codes
Main codes:
run_lpfc.m generates the necessary data files from monkey recording data.
run_lpfc.py generates the main figures.
supplementary.py generates the supplementary figures.

Other codes:
lpfc_data.py generates simulation data.
lpfc_model.py contains functions for simulations and analysis.
lpfc_para.py specifies model parameters (set overwrite to 1 if your simulation is interrupted).

Make sure that folders codes_matlab and Cai2014, as well as the above files are in the same folder.

### MATLAB
The MATLAB codes generate the necessary data files from monkey recording data for further analysis in Python. Modification of codes may be necessary for later version of MATLAB.

Data files needed in codes_matlab (contact person: X Cai): all X_dirtuningPar.mat (2662 files), X_psyphycell.mat (2662), X_tuning.mat (2662) and X_bin100.mat (2116, in profiles folder)

Matlab codes needed in Cai2014: arearead_DT.m, condReadout.m, readsession.m, sessionlistread.m, targetconvert.m, sessionlist_DT_PFC.m, degreediff.m, resultant.m, preferdir.m, hartigansdipsigniftest.m, hartigansdiptest.m

Files generated by Matlab: brainarea_all.txt (which LPFC region, v or d), dipp_all.txt, dipp_alld.txt, dipp_allv.txt, dipp_selected.txt, dipp_selectedd.txt, dipp_selectedv.txt (dip value and p-value for the peak difference distribution), peak_significant.txt (p-value of ranksum test: A0, B0, A400, B400), pval_CJ.txt (p-value of CJ encoding in windows 5 and 6, prior to target on), selectedID.txt (selected neurons), tuningallA0.txt, tuningallA400.txt, tuningallB0.txt, tuningallB400.txt, (8-direction tuning curves when A or B chosen in different windows), tuningA0,..., tuningA500, tuningB0,..., tuningB500 (same as above but for selected neurons)

You may need to modify the codes if you use a newer version of MATLAB.

### Python
The Python codes simulate the neural circuit models, analyze the data and generate the figures.

### License
This project is licensed under the MIT License.
