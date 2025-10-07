#!/bin/bash

# this script uses treatment region provided by Matt and frontal DMN network as control region

# FSL script for analyzing pre/post treatment resting state fMRI data
# This pipeline analyzes connectivity changes and activity differences

# ===== SETUP VARIABLES =====
# Edit these paths for your specific environment
STUDY_DIR="/export/neuroanaly-q3779/Projects/PilotUS/results/rsfmri3b/"
PRE_DATA_DIR="/export/neuroanaly-q3779/Projects/PilotUS/results/FreeSurfer_7_2_0__0_8_T1/CROSS/"
POST_DATA_DIR="${PRE_DATA_DIR}"
MASK_DIR="${STUDY_DIR}/masks"
OUTPUT_DIR="${STUDY_DIR}/results"
TR=0.735  # Repetition time in seconds



# Create output directory if it doesn't exist
mkdir -p ${STUDY_DIR}
mkdir -p ${OUTPUT_DIR}/individual
mkdir -p ${OUTPUT_DIR}/group


# Initialize headers for CSV files if they don't exist yet
if [ ! -f ${OUTPUT_DIR}/activity_changes.csv ]; then
  echo "Subject,TreatmentRegionChange,ControlRegionChange" > ${OUTPUT_DIR}/activity_changes.csv
fi

if [ ! -f ${OUTPUT_DIR}/roi_connectivity.csv ]; then
  echo "Subject,PreTreatmentCorrelation,PostTreatmentCorrelation" > ${OUTPUT_DIR}/roi_connectivity.csv
fi


# ===== PREPROCESSING =====
# Loop through each subject
for subject in 01-001 01-002 01-003 01-005 01-006 01-007 01-008 01-010 01-011 01-012 01-013; do
  echo "Processing subject: ${subject}"
  
  # Define subject directories
  PRE_SUBJECT_DIR="${PRE_DATA_DIR}/${subject}/BL/NativeSpace/rsfmri2/"
  POST_SUBJECT_DIR="${PRE_DATA_DIR}/${subject}/EOS/NativeSpace/rsfmri2/"
  
  
  # Output directories for this subject
  SUBJ_OUTPUT="${OUTPUT_DIR}/individual/${subject}"
  mkdir -p ${SUBJ_OUTPUT}
  
  # === PRE-TREATMENT DATA ===
  echo "Processing pre-treatment data"
  
  # Basic preprocessing (if not already done)
  # Motion correction
  mcflirt -in ${PRE_SUBJECT_DIR}/DC_rsFMRI_ND.nii.gz -out ${SUBJ_OUTPUT}/pre_mc -mats -plots
  
  # Brain extraction
  fslmaths ${SUBJ_OUTPUT}/pre_mc -Tmean ${SUBJ_OUTPUT}/pre_avg.nii.gz
  #bet ${SUBJ_OUTPUT}/pre_mc ${SUBJ_OUTPUT}/pre_brain -f 0.3 -R
  bet ${SUBJ_OUTPUT}/pre_avg ${SUBJ_OUTPUT}/pre_brain -f 0.3 -R -m 
  fslmaths ${SUBJ_OUTPUT}/pre_mc -mas ${SUBJ_OUTPUT}/pre_brain_mask.nii ${SUBJ_OUTPUT}/pre_brain
  
  # Spatial smoothing (6mm FWHM)
  fslmaths ${SUBJ_OUTPUT}/pre_brain.nii.gz -s 2.548 ${SUBJ_OUTPUT}/pre_smooth
  
  # High-pass temporal filtering (100s cutoff)
  fslmaths ${SUBJ_OUTPUT}/pre_smooth -bptf 50 -1 ${SUBJ_OUTPUT}/pre_filtered

  
  # creating mask
  fslmaths ${SUBJ_OUTPUT}/pre_brain -Tmean ${SUBJ_OUTPUT}/pre_avg.nii.gz
  nifty_linear_registration -f ${SUBJ_OUTPUT}/pre_avg.nii.gz -m   ${PRE_DATA_DIR}/${subject}/BL/FSNativeSpace/GM_Treatment_mask.nii.gz -t ${PRE_SUBJECT_DIR}/tt.tfm -r rigid -s slow -o   ${PRE_SUBJECT_DIR}/GM_Treatmentcontrol_mask_native.nii.gz -i NN
  mri_binarize --i  ${PRE_SUBJECT_DIR}/GM_Treatmentcontrol_mask_native.nii.gz   --match 1 --o ${SUBJ_OUTPUT}/pre_treatment_mask.nii.gz  2>&1 > /dev/null ;
  
  rm   ${PRE_SUBJECT_DIR}/aparcaseg.nii.gz  ${PRE_SUBJECT_DIR}/t1w_fmrinative.nii.gz
  nifty_linear_registration -f ${SUBJ_OUTPUT}/pre_avg.nii.gz -m   ${PRE_DATA_DIR}/${subject}/BL/mri/norm.mgz -t ${PRE_SUBJECT_DIR}/tt2.tfm -r rigid -s slow -o   ${PRE_SUBJECT_DIR}/t1w_fmrinative.nii.gz -i NN  2>&1 > /dev/null ;
  nifty_linear_registration -f ${SUBJ_OUTPUT}/pre_avg.nii.gz -m   ${PRE_DATA_DIR}/${subject}/BL/mri/aparc+aseg.mgz -t ${PRE_SUBJECT_DIR}/tt2.tfm -r rigid -s slow -o   ${PRE_SUBJECT_DIR}/aparcaseg.nii.gz -i NN  
  mri_binarize --i  ${PRE_SUBJECT_DIR}/aparcaseg.nii.gz  --match 1002 1026 2002  2026  --o ${SUBJ_OUTPUT}/pre_control_mask.nii.gz  2>&1 > /dev/null ; # caudal and rostral anterior congulate

  
  # === POST-TREATMENT DATA ===
  echo "Processing post-treatment data"
  
  # Repeat preprocessing for post-treatment data
  mcflirt -in ${POST_SUBJECT_DIR}/DC_rsFMRI_ND.nii.gz -out ${SUBJ_OUTPUT}/post_mc -mats -plots
#  bet ${SUBJ_OUTPUT}/post_mc ${SUBJ_OUTPUT}/post_brain -f 0.3 -R
  fslmaths ${SUBJ_OUTPUT}/post_mc -Tmean ${SUBJ_OUTPUT}/post_avg.nii.gz
  bet ${SUBJ_OUTPUT}/post_avg ${SUBJ_OUTPUT}/post_brain -f 0.3 -R -m 
  fslmaths ${SUBJ_OUTPUT}/post_mc -mas ${SUBJ_OUTPUT}/post_brain_mask.nii ${SUBJ_OUTPUT}/post_brain 
  
  fslmaths ${SUBJ_OUTPUT}/post_brain.nii.gz -s 2.548 ${SUBJ_OUTPUT}/post_smooth
  fslmaths ${SUBJ_OUTPUT}/post_smooth -bptf 50 -1 ${SUBJ_OUTPUT}/post_filtered


  # creating mask
  fslmaths ${SUBJ_OUTPUT}/post_brain -Tmean ${SUBJ_OUTPUT}/post_avg.nii.gz
  nifty_linear_registration -f ${SUBJ_OUTPUT}/post_avg.nii.gz -m   ${POST_DATA_DIR}/${subject}/EOS/FSNativeSpace/GM_Treatment_mask.nii.gz -t ${POST_SUBJECT_DIR}/tt.tfm -r rigid -s slow -o   ${POST_SUBJECT_DIR}/GM_Treatmentcontrol_mask_native.nii.gz  -i NN
  mri_binarize --i  ${POST_SUBJECT_DIR}/GM_Treatmentcontrol_mask_native.nii.gz   --match 1 --o ${SUBJ_OUTPUT}/post_treatment_mask.nii.gz  2>&1 > /dev/null ;

   rm   ${POST_SUBJECT_DIR}/aparcaseg.nii.gz  ${POST_SUBJECT_DIR}/t1w_fmrinative.nii.gz
  nifty_linear_registration -f ${SUBJ_OUTPUT}/post_avg.nii.gz -m   ${POST_DATA_DIR}/${subject}/EOS/mri/norm.mgz -t ${POST_SUBJECT_DIR}/tt2.tfm -r rigid -s slow -o   ${POST_SUBJECT_DIR}/t1w_fmrinative.nii.gz -i NN  2>&1 > /dev/null ;
  nifty_linear_registration -f ${SUBJ_OUTPUT}/post_avg.nii.gz -m   ${POST_DATA_DIR}/${subject}/EOS/mri/aparc+aseg.mgz -t ${POST_SUBJECT_DIR}/tt2.tfm -r rigid -s slow -o   ${POST_SUBJECT_DIR}/aparcaseg.nii.gz -i NN
  mri_binarize --i  ${POST_SUBJECT_DIR}/aparcaseg.nii.gz  --match 1002 1026 2002  2026 --o ${SUBJ_OUTPUT}/post_control_mask.nii.gz  2>&1 > /dev/null ; # caudal and rostral anterior congulate

  # ===== ACTIVITY ANALYSIS =====
  echo "Calculating activity metrics"

  
  # Extract mean time series from treatment and control regions (pre-treatment)
  fslmeants -i ${SUBJ_OUTPUT}/pre_filtered -o ${SUBJ_OUTPUT}/pre_treatment_timeseries.txt -m ${SUBJ_OUTPUT}/pre_treatment_mask.nii.gz
  fslmeants -i ${SUBJ_OUTPUT}/pre_filtered -o ${SUBJ_OUTPUT}/pre_control_timeseries.txt -m ${SUBJ_OUTPUT}/pre_control_mask.nii.gz
  
  # Extract mean time series from regions (post-treatment)
  fslmeants -i ${SUBJ_OUTPUT}/post_filtered -o ${SUBJ_OUTPUT}/post_treatment_timeseries.txt -m ${SUBJ_OUTPUT}/post_treatment_mask.nii.gz 
  fslmeants -i ${SUBJ_OUTPUT}/post_filtered -o ${SUBJ_OUTPUT}/post_control_timeseries.txt -m ${SUBJ_OUTPUT}/post_control_mask.nii.gz
  
  # Calculate temporal standard deviation (measure of activity fluctuation)
  fslmaths ${SUBJ_OUTPUT}/pre_filtered -Tstd ${SUBJ_OUTPUT}/pre_std
  fslmaths ${SUBJ_OUTPUT}/post_filtered -Tstd ${SUBJ_OUTPUT}/post_std
  
  # Extract mean standard deviation within masks
  fslstats ${SUBJ_OUTPUT}/pre_std -k ${SUBJ_OUTPUT}/pre_treatment_mask.nii.gz -m > ${SUBJ_OUTPUT}/pre_treatment_std.txt
  fslstats ${SUBJ_OUTPUT}/pre_std -k ${SUBJ_OUTPUT}/pre_control_mask.nii.gz -m > ${SUBJ_OUTPUT}/pre_control_std.txt
  fslstats ${SUBJ_OUTPUT}/post_std -k ${SUBJ_OUTPUT}/post_treatment_mask.nii.gz -m > ${SUBJ_OUTPUT}/post_treatment_std.txt
  fslstats ${SUBJ_OUTPUT}/post_std -k ${SUBJ_OUTPUT}/post_control_mask.nii.gz -m > ${SUBJ_OUTPUT}/post_control_std.txt
  
  # Calculate percentage change in activity (post - pre)/pre * 100
  pre_treatment_std=$(cat ${SUBJ_OUTPUT}/pre_treatment_std.txt)
  post_treatment_std=$(cat ${SUBJ_OUTPUT}/post_treatment_std.txt)
  pre_control_std=$(cat ${SUBJ_OUTPUT}/pre_control_std.txt)
  post_control_std=$(cat ${SUBJ_OUTPUT}/post_control_std.txt)
  
  treatment_pct_change=$(echo "scale=4; ($post_treatment_std - $pre_treatment_std) / $pre_treatment_std * 100" | bc)
  control_pct_change=$(echo "scale=4; ($post_control_std - $pre_control_std) / $pre_control_std * 100" | bc)
  
  echo "${subject},${treatment_pct_change},${control_pct_change}" >> ${OUTPUT_DIR}/activity_changes.csv
  
  # ===== FUNCTIONAL CONNECTIVITY ANALYSIS =====
  echo "Calculating functional connectivity metrics"
  
  # Seed-based connectivity (using treatment region as seed)
  fsl_glm -i ${SUBJ_OUTPUT}/pre_filtered -d ${SUBJ_OUTPUT}/pre_treatment_timeseries.txt -o ${SUBJ_OUTPUT}/pre_treatment_connectivity --demean
  fsl_glm -i ${SUBJ_OUTPUT}/post_filtered -d ${SUBJ_OUTPUT}/post_treatment_timeseries.txt -o ${SUBJ_OUTPUT}/post_treatment_connectivity --demean
  
  # Z-transform connectivity maps
  fslmaths ${SUBJ_OUTPUT}/pre_treatment_connectivity -add 1 -log -mul 0.5 ${SUBJ_OUTPUT}/pre_treatment_connectivity_z
  fslmaths ${SUBJ_OUTPUT}/post_treatment_connectivity -add 1 -log -mul 0.5 ${SUBJ_OUTPUT}/post_treatment_connectivity_z
  
  # Calculate connectivity difference (post - pre)
  fslmaths ${SUBJ_OUTPUT}/post_treatment_connectivity_z -sub ${SUBJ_OUTPUT}/pre_treatment_connectivity_z ${SUBJ_OUTPUT}/connectivity_diff
  
  # Get mean connectivity with control region
  fslmeants -i ${SUBJ_OUTPUT}/pre_treatment_connectivity_z -o ${SUBJ_OUTPUT}/pre_connectivity_treatment_to_control.txt -m ${SUBJ_OUTPUT}/pre_control_mask.nii.gz
  fslmeants -i ${SUBJ_OUTPUT}/post_treatment_connectivity_z -o ${SUBJ_OUTPUT}/post_connectivity_treatment_to_control.txt -m ${SUBJ_OUTPUT}/post_control_mask.nii.gz
  
  # ROI-to-ROI connectivity
  # Calculate Pearson correlation between treatment and control regions (pre)
  pre_treatment=$(cat ${SUBJ_OUTPUT}/pre_treatment_timeseries.txt)
  pre_control=$(cat ${SUBJ_OUTPUT}/pre_control_timeseries.txt)
  Rscript -e "cat(cor(scan('${SUBJ_OUTPUT}/pre_treatment_timeseries.txt'), scan('${SUBJ_OUTPUT}/pre_control_timeseries.txt')), '\n')" > ${SUBJ_OUTPUT}/pre_roi_to_roi_corr.txt
  
  # Calculate Pearson correlation between treatment and control regions (post)
  post_treatment=$(cat ${SUBJ_OUTPUT}/post_treatment_timeseries.txt)
  post_control=$(cat ${SUBJ_OUTPUT}/post_control_timeseries.txt)
  Rscript -e "cat(cor(scan('${SUBJ_OUTPUT}/post_treatment_timeseries.txt'), scan('${SUBJ_OUTPUT}/post_control_timeseries.txt')), '\n')" > ${SUBJ_OUTPUT}/post_roi_to_roi_corr.txt
  
  # Extract values for group analysis
  pre_roi_corr=$(cat ${SUBJ_OUTPUT}/pre_roi_to_roi_corr.txt)
  post_roi_corr=$(cat ${SUBJ_OUTPUT}/post_roi_to_roi_corr.txt)
  echo "${subject},${pre_roi_corr},${post_roi_corr}" >> ${OUTPUT_DIR}/roi_connectivity.csv
  
  # Register to standard space for group analysis (assumes you have T1 and transformations)
  # This step would require additional registration steps that depend on your specific setup
  # Example: flirt -in ${SUBJ_OUTPUT}/connectivity_diff -ref $FSLDIR/data/standard/MNI152_T1_2mm -out ${SUBJ_OUTPUT}/connectivity_diff_std -init ${PRE_SUBJECT_DIR}/reg/example_func2standard.mat -applyxfm
  
done 

