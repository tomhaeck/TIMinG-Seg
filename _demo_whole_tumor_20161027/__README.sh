#!/bin/bash

# Patient HG0001 from the MICCAI BRATS 2012-2013 training data set [1]

# STEP 1: BIAS FIELD CORRECTION WITH METHOD AND SOFTWARE OF YOUR CHOICE
# recommended:
#    ANTS N4BiasFieldCorrection [2]:
#    >> N4BiasFieldCorrection -i t1.nii -o t1_bias.nii -b [500,3]
#    >> N4BiasFieldCorrection -i t2.nii -o t2_bias.nii -b [500,3]
#    >> N4BiasFieldCorrection -i flair.nii -o flair_bias.nii -b [500,3]

# STEP 2: REGISTRATION OF ATLAS (MOVING) TO PATIENT IMAGE (FIXED) WITH METHOD AND SOFTWARE OF YOUR CHOICE
# recommended:
#    ANTS:
#    >> ANTS 3 -m CC[t1_bias.nii,../_atlas/T1.nii,1,4] -i 50x30x10 -o ab.nii
#    >> WarpImageMultiTransform 3 ../_atlas/WM.nii _WM_warped.nii -R t1_bias.nii abWarp.nii abAffine.txt
#    >> WarpImageMultiTransform 3 ../_atlas/GM.nii _GM_warped.nii -R t1_bias.nii abWarp.nii abAffine.txt
#    >> WarpImageMultiTransform 3 ../_atlas/CSF.nii _CSF_warped.nii -R t1_bias.nii abWarp.nii abAffine.txt

# STEP 3: LESION SEGMENTATION WITH TIMinG-Seg
TIMinG-Seg -i flair_bias.nii -i t2_bias.nii -p _WM_warped.nii -p _GM_warped.nii -p _CSF_warped.nii -o timingseg.nii --mu 1e1 --alpha 0.7 --mask brainmask.nii --hyp 0 --hyp 0 --hypamount 2.0  --hypamount 2.0 --inclusion -1 --maxIterEM 1





# [1] B. H. Menze, A. Jakab, S. Bauer, J. Kalpathy-Cramer, K. Farahani, J. Kirby, Y. Burren, N. Porz, J. Slotboom, and R. Wiest. “The multimodal brain tumor image segmentation benchmark (BRATS)”. In: IEEE TMI 34.10 (2015), pp. 1993–2024.
# [2] B. B. Avants, C. L. Epstein, M. Grossman, and J. C. Gee. “Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain”. In: Media 12.1 (2008), pp. 26–41.
# [2] http://stnava.github.io/ANTs/ (on October 2016)