'''
# For plotting single fiber Vm and Im at electrodes

elec_dist = simBunch.elec_dist*1000
internodal_len = subset[-1].params.internodal_len
num_nodes = subset[-1].params.num_nodes
idx_center_node = (num_nodes-1)/2   # also half-len
V1 = np.zeros(len(t_vec_py))
V2 = np.zeros(len(t_vec_py))
V3 = np.zeros(len(t_vec_py))

for i in range(num_nodes):
    node_offset = (i-idx_center_node)*internodal_len
    dist1 = np.sqrt((node_offset - e1_offset)**2 + elec_dist**2)
    dist2 = np.sqrt((node_offset - e2_offset)**2 + elec_dist**2)
    dist3 = np.sqrt((node_offset - e3_offset)**2 + elec_dist**2)
    node_row = center_row-idx_center_node+i # find correct entry in array
    for t in range(len(t_vec_py)):
        if Im_vec_py[node_row,t,0] is not np.nan:
            V1[t] = V1[t] + Im_vec_py[node_row, t, -1]/(4*pi*dist1)
            V2[t] = V2[t] + Im_vec_py[node_row, t, -1]/(4*pi*dist2)
            V3[t] = V3[t] + Im_vec_py[node_row, t, -1]/(4*pi*dist3)

plt.figure()
plt.plot(t_vec_py, V1, 'b-')
plt.plot(t_vec_py, V2, 'r-')
plt.plot(t_vec_py, V3, 'g-')
plt.show(block=False)
'''

'''
# For running a single fiber simulation
from PNS import *
import numpy as np
import pdb

# Initialize sweeney nerve bundles
initBunch = Bunch()
initBunch.scale_diam = 0.6
initBunch.num_nodes = 299
initBunch.num_seg = 1
initBunch.rho_a = 54.7      # [ohm-cm]
initBunch.rho_e = 500       # [ohm-cm]
initBunch.node_cm = 2.5     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0
initBunch.fiber_diam = 20.0

nerve1 = MyelinatedNerve(initBunch, sweeneyCh)
stim = h.Ipulse2(nerve1.axons[0](0.5))
stim.delay = 2   # [mS]
stim.dur = 0.1   # [mS]
stim.per = 1000.0/15 # period in [mS] ... 15 Hz
#stim.amp = 0.4 # [nA] for 3.48 um
stim.amp = 3.5 # [nA] for 20.0um
stim.num = 150   

tstop = 500 # [ms]
dt = 0.015  # [ms]
v_init = -80.0 # [mV]
celsius = 37.0 # [deg C]
num_timesteps = int(tstop/dt)+1
Vm_vec = []
Im_vec = []

for i in range(nerve1.params.num_nodes):
    Vm_vec.append(h.Vector(num_timesteps,0))
    Vm_vec[i].record(nerve1.axons[i](0.5)._ref_v)
    Im_vec.append(h.Vector(num_timesteps, 0))
    Im_vec[i].record(nerve1.axons[i](0.5)._ref_i_membrane)

t_vec = h.Vector()
t_vec.record(h._ref_t)

# start sim
h.dt = dt
h.tstop = tstop
h.celsius = celsius
h.finitialize(v_init)
while (h.t < h.tstop):
    h.fadvance()

# do plotting
t_vec_py = np.array(t_vec)
Vm_vec_py = np.array(Vm_vec)
I_vecPy = np.array(Im_vec)

nodes_name = [5, 25, 50, 75, 100, 150, 250]
f, axarr = plt.subplots(len(nodes_name), sharex=True)
for i in range(len(nodes_name)): 
    axarr[i].plot(t_vec_py, Vm_vec_py[i,:], 'b')
axarr[-1].set_xlabel('Time (mS)')
axarr[0].set_ylabel('Vm (mV)')
plt.show(block=True)
'''

'''
#----------------------
# Testing saving data
#----------------------
import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from bunch import Bunch

pop1 = [15, 2.5]    # alpha, beta motor
pop2 = [4, 1.5]     # pain fibers

np.random.seed(13)  # set seed so we get fixed results
bundle_diam = np.concatenate([np.random.normal(15, 2.5, 50), np.random.normal(4, 1.5, 50)])
base_freq = np.random.randint(10, 20, 100)
active_freq = np.random.randint(30, 60, 100)

initBunch = Bunch()
initBunch.scale_diam = 0.6
initBunch.num_nodes = 51    # number of nodes the longest fiber will have
initBunch.fiber_diam = 12   # default...will change
initBunch.num_seg = 1
initBunch.rho_a = 54.7      # [ohm-cm]
initBunch.rho_e = 500       # [ohm-cm]
initBunch.node_cm = 2.5     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0

max_internodal_len = initBunch.scale_len * bundle_diam[np.argmax(bundle_diam)]+0.5*initBunch.node_len
max_nerve_len = max_internodal_len * initBunch.num_nodes

# Get the number of nodes needed for each fiber
list_num_nodes = []
for i in np.arange(len(bundle_diam)):
    cur_internodal_len = initBunch.scale_len*bundle_diam[i] + 0.5*initBunch.node_len
    list_num_nodes.append(int(max_nerve_len/cur_internodal_len)//2*2+1)

# Try fiber 10 and 90 -- plot voltage over entire fiber every 100ms
# Expect to see Action potentials at each node to be slightly shifted from each other
# in time.
t_vec_py = np.load('./sim_data/t_vec_py_baseRate.npy')

Vm_vec_py_10 = np.load('./sim_data/Vm_vec_py_fiber10_baseRate.npy')
internodal_len_10 = initBunch.scale_len*bundle_diam[10] + 0.5*initBunch.node_len
n_nodes_10 = list_num_nodes[10]

Vm_vec_py_90 = np.load('./sim_data/Vm_vec_py_fiber90_baseRate.npy')
internodal_len_90 = initBunch.scale_len*bundle_diam[90] + 0.5*initBunch.node_len
n_nodes_90 = list_num_nodes[90]

# dt = 0.015ms...plot Vm at different nodes [0, 5, 15, 30] over all time
nodeIdx = [0, 5, 15, 30]

f, axarr = plt.subplots(len(nodeIdx), sharex=True)
for i in range(len(nodeIdx)):
    axarr[i].set_title('Node %d' % (nodeIdx[i]))
    axarr[i].set_ylabel('Vm (mV)')
    axarr[i].plot(t_vec_py, Vm_vec_py_10[nodeIdx[i], :])
axarr[-1].set_xlabel('Time (mS)')
axarr[0].set_title('Fiber 10, baseRate=%d Hz, diam=%.2f um' % (base_freq[10], bundle_diam[10]))
plt.show(block=False)

f, axarr = plt.subplots(len(nodeIdx), sharex=True)
for i in range(len(nodeIdx)):
    axarr[i].set_title('Node %d' % (nodeIdx[i]))
    axarr[i].set_ylabel('Vm (mV)')
    axarr[i].plot(t_vec_py, Vm_vec_py_90[nodeIdx[i], :])
axarr[-1].set_xlabel('Time (mS)')
axarr[0].set_title('Fiber 90, baseRate=%d Hz, diam=%.2f um' % (base_freq[90], bundle_diam[90]))
plt.show(block=False)
'''

#----------------------------------------------
# Test effects of electrode size on fiber's ENG
#----------------------------------------------
import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from bunch import Bunch

np.random.seed(13)  # set seed so we get fixed results
bundle_diam = np.concatenate([np.random.normal(15, 2.5, 50), np.random.normal(4, 1.5, 50)])
base_freq = np.random.randint(10, 20, 100)
active_freq = np.random.randint(30, 60, 100)

initBunch = Bunch()
initBunch.scale_diam = 0.6
initBunch.num_nodes = 51    # number of nodes the longest fiber will have
initBunch.fiber_diam = 12   # default...will change
initBunch.num_seg = 1
initBunch.rho_a = 54.7      # [ohm-cm]
initBunch.rho_e = 500       # [ohm-cm]
initBunch.node_cm = 2.5     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0

max_internodal_len = initBunch.scale_len * bundle_diam[np.argmax(bundle_diam)]+0.5*initBunch.node_len
max_nerve_len = max_internodal_len * initBunch.num_nodes

# Get the number of nodes needed for each fiber
list_num_nodes = []
for i in np.arange(len(bundle_diam)):
    cur_internodal_len = initBunch.scale_len*bundle_diam[i] + 0.5*initBunch.node_len
    list_num_nodes.append(int(max_nerve_len/cur_internodal_len)//2*2+1)

# Try fiber 10 and 90 -- 
# For each fiber finds the voltage measured at electrode spaced different
# distance apart and plot tripolar recordings.

fiberIdx = 90
rate = 'baseRate'
t_vec_py = np.load('./sim_data/t_vec_py_%s.npy' % (rate))

K = [4.0, 5.0, 6.0, 10.0, 14.0, 15.0]   # recording electrode spacing in [mm]
V1 = []
V2 = []
V3 = []

# Initialize electrode measurment
for kk in K:
    V1.append(np.zeros(len(t_vec_py)))
    V2.append(np.zeros(len(t_vec_py)))
    V3.append(np.zeros(len(t_vec_py)))

# Load fiber one-by-one and calculate contribution to electrode
print 'Nerve %d...' % (fiberIdx)
internodal_len = initBunch.scale_len * bundle_diam[fiberIdx] + 0.5*initBunch.node_len
num_nodes = list_num_nodes[fiberIdx]
idx_center_node = (num_nodes-1)/2   # also half-len

cur_Im_vec = np.load('./sim_data/Im_vec_py_fiber%d_%s.npy' % (fiberIdx, rate))
node_diam = initBunch.scale_diam * bundle_diam[fiberIdx]
area = pi * node_diam * initBunch.node_len
cur_Im_vec *= area

elec_dist = 1000
for kIdx in range(len(K)):
    kval = K[kIdx]
    e2_offset = 0              # electrode 2 aligned with center electrode
    e1_offset = -kval*(10**3)  # electrode 1 offset from center [um]
    e3_offset = kval*(10**3)   # electrode 3 offset from center [um]

    node_vec = np.arange(num_nodes)
    node_vec.shape = (num_nodes, 1) # conver to column vector
    node_offset = (node_vec - idx_center_node)*internodal_len
    dist1 = np.sqrt((node_offset - e1_offset)**2 + elec_dist**2)
    dist2 = np.sqrt((node_offset - e2_offset)**2 + elec_dist**2)
    dist3 = np.sqrt((node_offset - e3_offset)**2 + elec_dist**2)

    V1[kIdx] += np.sum(cur_Im_vec/(4*pi*dist1), 0)
    V2[kIdx] += np.sum(cur_Im_vec/(4*pi*dist2), 0)
    V3[kIdx] += np.sum(cur_Im_vec/(4*pi*dist3), 0)

'''
    # Vectorized code is almost 100 times faster..
    for i in range(num_nodes):
        node_offset = (i-idx_center_node)*internodal_len
        dist1 = np.sqrt((node_offset - e1_offset)**2 + elec_dist**2)
        dist2 = np.sqrt((node_offset - e2_offset)**2 + elec_dist**2)
        dist3 = np.sqrt((node_offset - e3_offset)**2 + elec_dist**2)

        for t in range(len(t_vec_py)):
            V1[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist1)
            V2[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist2)
            V3[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist3)
'''
            
# Plot measurement for different K-placements
Vpp = np.zeros(len(K))
f, axarr = plt.subplots(len(K), sharex=True)
for kIdx in range(len(K)):
    V_ref = (V1[kIdx] + V3[kIdx])/2
    Vsig = V2[kIdx] - V_ref
    Vpp[kIdx] = max(Vsig)-min(Vsig)
    axarr[kIdx].plot(t_vec_py, Vsig)
    axarr[kIdx].set_title('K=%d mm' % (K[kIdx]))
    axarr[kIdx].set_ylabel('V (mV)')

axarr[0].set_title('%s, fiber size=%.2f um, K=%d mm' % (rate, bundle_diam[fiberIdx], K[0]))
axarr[-1].set_xlabel('Time (mS)')
plt.show(block=False)

# Plot Vpp as a function of K
plt.figure()
plt.plot(K, Vpp, 'b-*')
plt.xlabel('K (mm)')
plt.ylabel('Vpp (mV)')
plt.show(block=False)
