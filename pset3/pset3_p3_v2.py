from PNS import *
import numpy as np

'''
Not enough memory...
Do simulation fiber by fiber independently, then save the data
'''
pop1 = [15, 2.5]    # alpha, beta motor
pop2 = [4, 1.5]     # pain fibers

np.random.seed(13)  # set seed so we get fixed results
bundle_diam = np.concatenate([np.random.normal(15, 2.5, 50), np.random.normal(4, 1.5, 50)])
base_freq = np.random.randint(10, 20, 100)

# Initialize sweeney nerve bundles
initBunch = Bunch()
initBunch.scale_diam = 0.6
initBunch.num_nodes = 51
initBunch.num_seg = 1
initBunch.rho_a = 54.7      # [ohm-cm]
initBunch.rho_e = 500       # [ohm-cm]
initBunch.node_cm = 2.5     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0

for i in range(len(bundle_diam)):
    max_internodal_len = initB
