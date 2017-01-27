import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from neuron import h, gui

# Make sure all parameters need to be initialized as doubles for accurate math!!

# Independent params
fiber_diam = 1.0
num_nodes = 101
my_nseg = 3

# Derived params
node_diam = 0.7*fiber_diam           # node diameter [um]
node_diam_cm = node_diam/10**4       # node diameter [cm]
node_length = 1.5                    # node length [um]
node_length_cm = node_length/10**4   # node length [cm]
myelin_length = 100.0                  # internodal myelin length [um]
myelin_length_cm = myelin_length/10**4 
elec_dist = 1.0                        # shortest dis from electrode to a node [mm]
elec_dist_um = elec_dist*10**3       # convert to neuron default length unit [um]

# For calculating extracellular V due to electrode [um]
internodal_length = myelin_length + 0.5*node_length

# Electrical properties
node_cm = 1.0					# specific membrane capacitance [uF/cm**2]
node_rm = 2000.0              # specific membrane resistance [ohm-cm**2]
rho_a = 100.0					# intracellular resistivity [ohm-cm]
rho_a_adj = rho_a*((node_length_cm+myelin_length_cm)/node_length_cm)    # adjusted rho_a so Neuron can get the same internodal resistance [ohm-cm]
rho_e = 500.0                 # extracellular resistivity [ohm-cm]
#rho_e_adj = ??

# Segment electrical properties -- for part 1 of problem 1
my_Ra = 4*rho_a*(myelin_length_cm+node_length_cm)/(pi*node_diam_cm**2)  # internodal resistance [ohm]
my_Cm = node_cm*(pi*node_length_cm)*node_diam_cm   # nodal membrane cap [uF]
my_Rm = node_rm/(pi*node_length_cm*node_diam_cm)   # nodal membrane res [ohm]

# Stimulus parameters
mydel = 5.0					# start at t=5ms [ms]
myamp = -5  				# amplitude [nA]. Cathodic current -- negative
mydur = 0.1					# duration, aka pulsewidth [ms] -- 0.1ms=100us

#------------ Params used in Neuron environment ------
# Temporal parameters
h.dt = 0.005					# [ms]
h.tstop = 10.0					# [ms]
num_timesteps = int(h.tstop/h.dt) + 1   # derived for local use

# Other parameters
h.v_init = -70.0				# [mV]
h.celsius = 37.0				# [deg C]

#--------------- Model Initialization ---------------
# make axons
axons = []
for i in range(num_nodes):
    axons.append(h.Section('axon'))
    axons[i].push()
    axons[i].nseg = my_nseg
    axons[i].diam = node_diam
    axons[i].L = node_length
    axons[i].Ra = rho_a_adj   # adjusted axoplasmic resistivity in ohm*cm
    axons[i].cm = node_cm 

    axons[i].insert('pas')
    axons[i].g_pas = 1.0/node_rm

    axons[i].insert('extracellular')
    h.pop_section()

# connect the axon sections
for i in range(num_nodes-2):
        axons[i+1].connect(axons[i], 1, 0)

#-------------- Instrumentation ----------------
dummy_elec = h.Section('dummy_elec')
dummy_stim = h.IClamp(dummy_elec(0.5))
dummy_stim.delay = mydel
dummy_stim.amp = myamp
dummy_stim.dur = mydur

def calcDist(curIdx, elecIdx):
    return np.sqrt((np.abs(curIdx-elecIdx)*internodal_length)**2 + elec_dist_um**2)

# calc distance from electrode to each axon
elec_node_idx = num_nodes/2
r_vec = [calcDist(i, elec_node_idx) for i in range(num_nodes)]

# record Vm(t) at all nodes
Vm_vec = []
for i in range(num_nodes):
    Vm_vec.append(h.Vector(num_timesteps,0))
    Vm_vec[i].record(axons[i](0.5)._ref_v)

# time vector
t_vec = h.Vector()
t_vec.record(h._ref_t)

#-------------- Simulation control ------------
def run_sim():
    h.finitialize(h.v_init)
    while(h.t < h.tstop):
        for i in range(num_nodes):
            # apply extracellular potential -- 10^3*(uV) = [mV]
            axons[i].e_extracellular = (10**3)*rho_e*dummy_stim.i/(4*pi*r_vec[i])
        h.fadvance()
    Vm_vec_py = np.array(Vm_vec)
    return Vm_vec_py

Vm_vec_py = run_sim()
while Vm_vec_py.max() < h.v_init+20:
    dummy_stim.amp = dummy_stim.amp-0.5
    Vm_vec_py = run_sim()
print "At %.2fnA extra-cellular stimulation, max change in Vm=%.2fmV" % (dummy_stim.amp, Vm_vec_py.max()-h.v_init)

#------------- Plot Data --------------
def plot_data(Vm_vec_py):
    #Vm_vec_py = np.array(Vm_vec)

    # Plot Vm(t, x=25)
    plt.figure(figsize=(8,4))
    plt.plot(t_vec, Vm_vec_py[elec_node_idx,:])
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.ylim((-80,-30))
    plt.show(block=False)
    plt.title("Vm(t, x=%d)" % (elec_node_idx))

    # Plot Vm(t=5.5ms, x)
    plt.figure(figsize=(8,4))
    tStart = 4.0
    tEnd = 7.0
    for i in np.arange(tStart, tEnd, 0.1): 
        v_alpha = (i-tStart)/(tEnd-tStart)
        plt.plot(np.arange(num_nodes), Vm_vec_py[:, int(i/h.dt)], 'r-', alpha=v_alpha)
        plt.xlabel('node')
        plt.ylabel('mV')
        plt.ylim((-80,-30))
    plt.title("Vm(t=, x)")
    plt.show(block=False)
