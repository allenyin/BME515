"""
BME515 Pset 1
"""
import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from neuron import h, gui

def um_to_cm(x):
    return x / (10**4)

def cm_to_um(x):
    return x * (10**4)

def calcParams(D):
    """ 
    Return cell properties as a function of D, in um
    """
    # Geometrical properties
    node_diam = 0.7*D       # node diameter [um]
    node_length = 1.5       # node length [um]
    myelin_length = 100*D   # internodal length [um]

    # Electrical properties
    node_cm = 1             # specific mem capacitance [uF/cm^2]
    node_rm = 2000          # specific mem resistance [ohm*cm^2]
    rho_a = 100             # intracellular resistivity [ohm*cm]

    # Calculate segment properties
    # First lengths to cm
    node_diam_cm = um_to_cm(node_diam)
    node_length_cm = um_to_cm(node_length)
    myelin_length_cm = um_to_cm(node_length)
    
    Ra = 4*rho_a*(myelin_length_cm+node_length_cm)/(pi * node_diam_cm**2)  # internodal resistance [ohm]...is this resistivity right?...
    Cm = node_cm*(pi*node_length_cm)*node_diam_cm   # nodal mem cap [uF]
    Rm = node_rm/(pi*node_length_cm*node_diam_cm)   # nodal mem res [ohm]

    rho_a_adj = rho_a*((node_length_cm+myelin_length_cm)/node_length_cm)  # adjusted rho_a for Neuron, [ohm*cm]

    geometry_params = {}
    geometry_params['node_diam'] = node_diam
    geometry_params['node_length'] = node_length
    geometry_params['myelin_length'] = myelin_length
    geometry_params['node_cm'] = node_cm
    geometry_params['node_rm'] = node_rm
    geometry_params['rho_a'] = rho_a
    geometry_params['rho_a_adj'] = rho_a_adj
    geometry_params['node_diam_cm'] = node_diam_cm
    geometry_params['node_length_cm'] = node_length_cm
    geometry_params['myelin_length_cm'] = myelin_length_cm
    geometry_params['Ra'] = Ra
    geometry_params['Cm'] = Cm
    geometry_params['Rm'] = Rm
    return geometry_params

def p1_1():
    Ra = []
    Cm = []
    Rm = []
    D = []
    for i in np.arange(1,21,2):
        params = calcParams(i)
        D.append(i)
        Ra.append(params['Ra'])
        Cm.append(params['Cm'])
        Rm.append(params['Rm'])

    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].plot(D, Ra)
    axarr[0].set_title('Internodal Axoplasmic resistance')
    axarr[0].set_ylabel('Ohm')

    axarr[1].plot(D, Cm)
    axarr[1].set_title('Nodal Membrane capacitance')
    axarr[1].set_ylabel('uF')

    axarr[2].plot(D, Rm)
    axarr[2].set_title('Nodal Membrane resistance')
    axarr[2].set_ylabel('Ohm')
    axarr[2].set_xlabel('D (um)')

    plt.show(block=False)
    plt.savefig('p1_1.png', bbox_inches='tight')

#-----------------------------------------------------------
# Neuron simulation start
#-----------------------------------------------------------
def init_sim():
    h.finitialize(-70)  # init voltage state [mV]
    h.celsius = 37            # set temp [deg C]
    h.tstop = 300             # simulation time [ms]

axons = []
def initialize_axons(D, nNodes, nseg):
    params = calcParams(D)
    for i in range(nNodes-1):
        axons.append(h.Section('axon'))
        axons[i].push()
        axons[i].nseg = nseg
        axons[i].diam = params['node_diam']
        axons[i].L = params['node_length']
        axons[i].Ra = params['rho_a_adj']   # adjusted axoplasmic resistivity in ohm*cm
        axons[i].cm = params['node_cm'] 

        axons[i].insert('pas')
        axons[i].g_pas = 1.0/params['node_rm']
        h.pop_section()

    # now connect the nodes
    for i in range(nNodes-2):
        axons[i+1].connect(axons[i], 1, 0)

#----------------------------------------------------------
# Instrumentation
#----------------------------------------------------------
dummy_elec = []
dummy_stim = []
def make_elec():
    # Create a section called dummy_elec and attach dummy_stim to it
    dummy_elec = h.Section('dummy_elec')
    dummy_stim = h.IClamp(dummy_elec(0.5))
    dummy_stim.delay = 5
    dummy_stim.amp = -6500
    dummy_stim.dur = 0.1
    return (dummy_elec, dummy_stim)
dummy_elec, dummy_stim = make_elec()

