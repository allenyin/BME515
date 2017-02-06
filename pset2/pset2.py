import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from neuron import h, gui
from bunch import Bunch
np.seterr(all='raise')

""" 
Parameters to vary script wise:
    fiber_diam, num_nodes, my_nseg, elec_dist, myamp, mydur

Sections/variables to keep track:
    axons, dummy_elec, dummy_stim, Vm_vec, t_vec

Use Sweeney channels. Modified in create_axons(), edit_one_axons()
"""

def calcParams(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, mydur=0.1):
    """
    Contains all numerical parameters, either local or global
    in Neuron's environment

    Inputs:
        fiber_diam: D from pset, in [um]
        num_nodes: number of axon nodes
        my_nseg: number of segments per section
        elec_dist: perp distance from electrode to axon [mm]
        myamp: stimulation current, cathodic=negative. [nA].

    Derived gemoetrical params according to Sweeney1987, requires changing:
    : node_diam (0.7*fiber_diam -> 0.6*fiber_diam)
    : v_init (-70mV -> -80mV)
    : rho_a (-100 Ohm-cm -> 54.7 Ohm-cm)
    : node_cm (1.0 uF/cm^2 -> 2.5uF/cm^2)
    """

    # Derived geometrical params
    node_diam = 0.6*fiber_diam           # node diameter [um]
    node_diam_cm = node_diam/10**4       # node diameter [cm]
    node_length = 1.5                    # node length [um]
    node_length_cm = node_length/10**4   # node length [cm]
    myelin_length = 100.0*fiber_diam     # internodal myelin length [um]
    myelin_length_cm = myelin_length/10**4 
    elec_dist_um = elec_dist*10**3       # convert [mm] to [um]

    # For calculating extracellular V due to electrode [um]
    internodal_length = myelin_length + 0.5*node_length

    # Electrical properties
    node_cm = 2.5					# specific membrane capacitance [uF/cm**2]
    node_rm = 2000.0                                    # specific membrane resistance for passive channel [ohm-cm**2]
    rho_a =  54.7					# intracellular resistivity [ohm-cm]
    rho_a_adj = rho_a*((node_length_cm+myelin_length_cm)/node_length_cm)    # adjusted rho_a so Neuron can get the same internodal resistance [ohm-cm]
    rho_e = 500.0                                       # extracellular resistivity [ohm-cm]

    # Segment electrical properties -- for part 1 of problem 1
    my_Ra = 4*rho_a*(myelin_length_cm+node_length_cm)/(pi*node_diam_cm**2)  # internodal resistance [ohm]
    my_Cm = node_cm*(pi*node_length_cm)*node_diam_cm   # nodal membrane cap [uF]
    my_Rm = node_rm/(pi*node_length_cm*node_diam_cm)   # nodal membrane res [ohm]

    # Stimulus parameters
    mydel = 5.0					# start at t=5ms [ms]
    #mydur = 0.1				# duration, aka pulsewidth [ms] -- 0.1ms=100us

    elec_node_idx = int(num_nodes/2)

    # simulation-control parameters
    dt = 0.005                          # [ms]
    tstop = 10.0                        # [ms]
    num_timesteps = int(h.tstop/h.dt)+1 # unitless
    v_init = -80.0                     # [mV]
    celsius = 37.0                      # [deg C]
    return Bunch(locals())

def editAxons(axons, params):
    if axons is None:
    # Create axons if empty
        return createAxons(params)
    
    # edit the existing sections
    for i in range(params.num_nodes):
        axons[i].push()
        axons[i] = editOneAxon(axons[i], params)
        h.pop_section()
    
    # add more sections if we don't have enough
    if len(axons) < params.num_nodes:
        extra_axons = createAxons(params)
        extra_axons[0].connect(axons[-1],1,0)   # connect them together
        axons.append(extra_axons)
    return axons
    
def editOneAxon(axon, params):
# Edit just one axon section
    axon.nseg = params.my_nseg
    axon.diam = params.node_diam
    axon.L = params.node_length
    axon.cm = params.node_cm
    axon.Ra = params.rho_a_adj

    # no more passive channels
    #axon.g_pas = 1.0/params.node_rm

    # Sweeney channels
    axon.ena = 35.64    # [mV]
    #axon.gnabar
    #axon.el
    #axon.gl
    return axon

def createAxons(params):
    axons = []
    for i in range(params.num_nodes):
        axons.append(h.Section('axon'))
        axons[i].push()
        axons[i].insert('sweeney')
        axons[i].insert('extracellular')
        axons[i] = editOneAxon(axons[i], params)
        h.pop_section()
    for i in range(params.num_nodes-2):
        axons[i+1].connect(axons[i],1,0)
    return axons

def editStim(dummy_elec, dummy_stim, params):
    if (dummy_elec is None) or (dummy_stim is None):
        dummy_elec = h.Section('dummy_elec')
        dummy_stim = h.IClamp(dummy_elec(0.5))
    # Edit current params
    dummy_stim.delay = params.mydel
    dummy_stim.amp = params.myamp
    dummy_stim.dur = params.mydur
    return (dummy_elec, dummy_stim)

def calcDist(curIdx, params):
    elecIdx = params.elec_node_idx
    internodal_length = params.internodal_length
    elec_dist_um = params.elec_dist_um
    return np.sqrt((np.abs(curIdx-elecIdx)*internodal_length)**2 + elec_dist_um**2) # [um]

def setup_vectors(axons, Vm_vec, t_vec, params):
    if Vm_vec is None:
        Vm_vec = create_Vm_vec(axons, params)
    elif Vm_vec[0].size() != params.num_timesteps:
        Vm_vec = [v.resize(params.num_timesteps) for v in Vm_vec]
    if t_vec is None:
        t_vec = h.Vector()
        t_vec.record(h._ref_t)
    return (Vm_vec, t_vec)

def create_Vm_vec(axons, params):
    Vm_vec = []
    for i in range(params.num_nodes):
        Vm_vec.append(h.Vector(params.num_timesteps,0))
        Vm_vec[i].record(axons[i](0.5)._ref_v)
    return Vm_vec

def run_sim(axons, dummy_stim, r_vec, params):
    h.dt = params.dt
    h.tstop = params.tstop
    h.celsius = params.celsius
    h.finitialize(params.v_init)
    while (h.t < h.tstop):
        for i in range(params.num_nodes):
            """
            apply extracellular potential -- 
            If dummy_stim.i is in mA, then
             I/(4*pi*sigma*r) = (I*rho_e)/(4*pi*r) = ([mA][ohm-cm])/([um]) -> multiply by 10^4 to get mV

            If dummy_stim.i is in nA, then
             I/(4*pi*sigma*r) = (I*rho_e)/(4*pi*4) = ([nA][ohm-cm])/([um]) -> multiply by (1/10^6)/(1/10^4)=10^-2 to get mV
            """
            axons[i].e_extracellular = (10**-2)*params.rho_e*dummy_stim.i/(4*pi*r_vec[i])
        h.fadvance()

def run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, mydur=0.1, axons=None, dummy_elec=None, dummy_stim=None, Vm_vec=None, t_vec=None):
    
    # generate all the parameters
    params = calcParams(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, mydur)    
    axons = editAxons(axons, params)
    (dummy_elec, dummy_stim) = editStim(dummy_elec, dummy_stim, params)
    r_vec = [calcDist(i, params) for i in range(params.num_nodes)]
    (Vm_vec, t_vec) = setup_vectors(axons, Vm_vec, t_vec, params)
    run_sim(axons, dummy_stim, r_vec, params)
    return (axons, dummy_elec, dummy_stim, Vm_vec, t_vec)

# Start doing the work
#plt.close('all')

# Initialize Neuron objects references
axons = None
dummy_elec = None
dummy_stim = None
Vm_vec = None
t_vec = None

#-------------------- Hwk Part --------------------------------
"""
Experiment...what does the action potential look like?

Use the criteria if >90% of the nodes reach potential > +8mV),
then we have achieved propagating action potential for cathodic stimulation.
"""
"""
fiber_diam = 12 # 12um
num_nodes = 51
my_nseg = 1
elec_dist = 1   # 1mm
myamp = -0.115e-3/1e-9
mydur = 0.20     # ms
(axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
        run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, mydur,\
                   axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
Vm_vec_py = np.array(Vm_vec)

# Plot time course along fiber during stimulation
plt.figure()
for i in np.arange(1000, 2000, 10):
    plt.plot(np.arange(1,52), Vm_vec_py[:,i], 'r-')
plt.show(block=False)
"""


def AP_status(Vm_vec_py):
    """
    Action potential determination:
    90% of the nodes' maximum Vm > +8mV
    The global maximum needs to be less than 13.8mV
    """
    try:
        AP = ((sum(np.max(Vm_vec_py, 1) >= 8)/(Vm_vec_py.shape[0]*1.0)) >= 0.8)
    except:
        import pdb; pdb.set_trace()
        print "here"
    tooHigh = (np.max(Vm_vec_py) - 13 > 1)

    if not AP:
        return -1
    elif AP and not tooHigh:
        return 0
    else:
        return 1

def p1_1(axons, dummy_elec, dummy_stim, Vm_vec, t_vec):
    """
    Fix electrode at 1mm, D=12um. 
    Increase duration from 10us to 10ms in steps of
    1, 2, 5.

    Find cathodic current required to achieve action potential
    """
    threshCur = []
    mydur = np.array([10e-6, 20e-6, 50e-6, \
                      100e-6, 200e-6, 500e-6, \
                      1e-3, 2e-3, 5e-3, 10e-3])     # in seconds
    mydur = mydur/1e-3  # convert to ms
    fiber_diam = 12     # 12um diameter
    num_nodes = 51      # 51 nodes
    my_nseg = 1         # 1 segment
    elec_dist = 1       # 1mm electrode dist

    lo_amp = -0.05e-3/1e-9  # low cathodic current limit [nA]
    hi_amp = -3e-3/1e-9     # high cathodic current limit [nA]

    for dur in mydur:
        if len(threshCur) > 0:
            hi_amp = threshCur[-1]

        myamp = np.mean((lo_amp, hi_amp))
        (axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
        run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, dur, \
                   axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
        Vm_vec_py = np.array(Vm_vec)
        status = AP_status(Vm_vec_py)

        while status != 0:
            if status == 1: # current too high
                hi_amp = myamp
                myamp = np.mean((lo_amp, myamp))
            elif status == -1: # current too low
                lo_amp = myamp
                myamp = np.mean((myamp, hi_amp))

            (axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
            run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, dur, \
                        axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
            Vm_vec_py = np.array(Vm_vec)
            status = AP_status(Vm_vec_py)

        print "D=%.2fum, r=%.2fmm, dur=%.2fms: %.2fnA or %.2fmA see action potential" % \
            (fiber_diam, elec_dist, dur, dummy_stim.amp, dummy_stim.amp/10**6)
        threshCur.append(dummy_stim.amp)

    # plot threshold current
    #import pdb; pdb.set_trace()
    plt.figure()
    plt.semilogx(mydur, np.abs(threshCur)/10**6, '*-')
    plt.xlabel('Duration (ms)')
    plt.ylabel('Threshold current magnitude (mA)')
    plt.title("Active Sweeney Strength-Duration")
    plt.show(block=False)
    plt.savefig('sweeney_strengthDuration.png', bbox_inches='tight')
    return (axons, dummy_elec, dummy_stim, Vm_vec, t_vec)

(axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
        p1_1(axons, dummy_elec, dummy_stim, Vm_vec, t_vec)



