## Multi Source Self Calibration
######Multi-Source Self-Calibration (MSSC) is a direction-dependent, calibration technique which provides an additional step to standard phase referencing. MSSC uses multiple faint sources detected within the primary beam and combines them together. The combined response of many sources across the field-of-view is more than sufficient to allow phase corrections to be derived. Each source have their CLEAN models divided into the visibilities which results in multiple point sources. These are stacked in the uv plane to increase the S/N, permitting self-calibration to become feasible. It is worth noting that this process only applies to wide-field VLBI data sets that detect and image multiple sources within one epoch.  Recent improvements in the capabilities of VLBI correlators is ensuring that wide-field VLBI is a reality and as a result there will be an increased number of experiments which can utilise MSSC. If MSSC is used, please reference Radcliffe et al. (submitted) 

######Pre-requisites: Parseltongue, python 2.7 & AIPS 
######Usage: Sources which have been split and averaged to an individual file per source need to be in the same directory as multi_source_self_cal_v?.?.py and MSSC_input.txt. After changing the input file, MSSC can be run on using the following commands ParselTongue multi_source_self_cal_v?.?.py MSSC_input.txt and it should run. 


###### REMEMBER TO CHECK THE AIPS DIRECTORY BEFORE RUNNING MSSC! IT WILL DELETE FILES ON THE DIRECTORY SPECIFIED IN THE INPUTS
