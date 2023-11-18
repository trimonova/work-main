#his_filter=['DYKE','FOLD','DYKE'] # FAULT SHEAR-ZONE FOLD TILT UNCONFORMITY DYKE PLUG default as loaded is to allow model sequences
his_filter=['FOLD','FOLD','SHEAR-ZONE']
display_number=1 # number of randomly selected models to display
from Noddyverse import display_models
import numpy as np

#%matplotlib inline
model = display_models(his_filter,display_number)
perm_matrix = np.ones((N_r_full, M_fi_full))