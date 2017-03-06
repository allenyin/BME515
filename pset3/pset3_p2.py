from PNS import *
import numpy as np
import pdb

'''
2 populations of fibers:
    myelinated motor (alpha, beta): mean=15um, sd=2.5um. So 2 s.d. = 5um
    pain fibers (type-III): mean=4um, sd=1.5um. So 2 s.d. = 3um.

    Randomly select which population a fiber belongs to, then draw from
    the corresponding population
'''
pop1 = [15, 2.5]    # alpha, beta motor
pop2 = [4, 1.5]     # pain fibers

np.random.seed(13)  # set seed so we get fixed results
bundle_diam = np.concatenate([np.random.normal(15, 2.5, 50), np.random.normal(4, 1.5, 50)])

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

nerveBundle = []
for i in np.arange(100):
    initBunch.fiber_diam = bundle_diam[i]   # [um]
    nerveBundle.append(MyelinatedNerve(initBunch, sweeneyCh))

# Initialize stimulating electrode
myamp = -0.15e-3/1e-9   # 0.15mA cathodic current for 0.2ms duration for 12um fiber
mydur = 0.4 # ms
mydel = 2
elec = DummyElectrode(mydel, myamp, mydur)

# Initialize bundle simulation




