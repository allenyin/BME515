import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from neuron import h, gui
from bunch import Bunch

""" 
Parameters to vary script wise:
    fiber_diam, num_nodes, my_nseg, elec_dist

Sections/variables to keep track:
    axons, dummy_elec, dummy_stim, Vm_vec, t_vec
"""

def calcParams(fiber_diam, num_nodes, my_nseg, elec_dist, myamp):
    """
    Contains all numerical parameters, either local or global
    in Neuron's environment

    Inputs:
        fiber_diam: D from pset, in [um]
        num_nodes: number of axon nodes
        my_nseg: number of segments per section
        elec_dist: perp distance from electrode to axon [mm]
        myamp: stimulation current, cathodic=negative. [nA].
    """

    # Derived geometrical params
    node_diam = 0.7*fiber_diam           # node diameter [um]
    node_diam_cm = node_diam/10**4       # node diameter [cm]
    node_length = 1.5                    # node length [um]
    node_length_cm = node_length/10**4   # node length [cm]
    myelin_length = 100.0*fiber_diam     # internodal myelin length [um]
    myelin_length_cm = myelin_length/10**4 
    elec_dist_um = elec_dist*10**3       # convert [mm] to [um]

    # For calculating extracellular V due to electrode [um]
    internodal_length = myelin_length + 0.5*node_length

    # Electrical properties
    node_cm = 1.0					# specific membrane capacitance [uF/cm**2]
    node_rm = 2000.0                                    # specific membrane resistance [ohm-cm**2]
    rho_a = 100.0					# intracellular resistivity [ohm-cm]
    rho_a_adj = rho_a*((node_length_cm+myelin_length_cm)/node_length_cm)    # adjusted rho_a so Neuron can get the same internodal resistance [ohm-cm]
    rho_e = 500.0                                       # extracellular resistivity [ohm-cm]

    # Segment electrical properties -- for part 1 of problem 1
    my_Ra = 4*rho_a*(myelin_length_cm+node_length_cm)/(pi*node_diam_cm**2)  # internodal resistance [ohm]
    my_Cm = node_cm*(pi*node_length_cm)*node_diam_cm   # nodal membrane cap [uF]
    my_Rm = node_rm/(pi*node_length_cm*node_diam_cm)   # nodal membrane res [ohm]

    # Stimulus parameters
    mydel = 5.0					# start at t=5ms [ms]
    mydur = 0.1					# duration, aka pulsewidth [ms] -- 0.1ms=100us
    elec_node_idx = int(num_nodes/2)

    # simulation-control parameters
    dt = 0.005                          # [ms]
    tstop = 10.0                        # [ms]
    num_timesteps = int(h.tstop/h.dt)+1 # unitless
    v_init = -70.0                      # [mV]
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
    axon.Ra = params.rho_a_adj
    axon.cm = params.node_cm
    axon.g_pas = 1.0/params.node_rm
    return axon

def createAxons(params):
    axons = []
    for i in range(params.num_nodes):
        axons.append(h.Section('axon'))
        axons[i].push()
        axons[i].insert('pas')
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

def run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, axons=None, dummy_elec=None, dummy_stim=None, Vm_vec=None, t_vec=None):
    
    # generate all the parameters
    params = calcParams(fiber_diam, num_nodes, my_nseg, elec_dist, myamp)    
    axons = editAxons(axons, params)
    (dummy_elec, dummy_stim) = editStim(dummy_elec, dummy_stim, params)
    r_vec = [calcDist(i, params) for i in range(params.num_nodes)]
    (Vm_vec, t_vec) = setup_vectors(axons, Vm_vec, t_vec, params)
    run_sim(axons, dummy_stim, r_vec, params)
    return (axons, dummy_elec, dummy_stim, Vm_vec, t_vec)

# Start doing the work
plt.close('all')

# Initialize Neuron objects references
axons = None
dummy_elec = None
dummy_stim = None
Vm_vec = None
t_vec = None

#params = calcParams(12, 51, 1, 1, 0.1)
#axons = editAxons(axons, params)

def p1_2(axons, dummy_elec, dummy_stim, Vm_vec, t_vec):
    """
    Fix electrode at 1mm from axon, vary D from 1um to 21um in 2um increment.
    Increase cathodic current until Vm=20mV. Record and plot the D vs. current
    relationship
    """
    threshCur = []
    D = np.arange(1.0, 22.0, 2.0)
    for d in D:
        fiber_diam = d    # starting d, in [um]
        num_nodes = 51    # must be int
        my_nseg = 1       # must be int
        elec_dist = 1.0   # starting 1mm
                
        lo_amp = -10**2     # low cathodic current limit [nA]
        hi_amp = -10**8     # high cathodic current  limit [nA]
        
        # starting in the middle -- binary search
        myamp = np.mean((lo_amp, hi_amp))      
        
        (axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
                run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, \
                           axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
        Vm_vec_py = np.array(Vm_vec)
        targetV = h.v_init+20
        
        while np.abs(Vm_vec_py.max() - targetV) > 0.1:
            if Vm_vec_py.max() > targetV:
                hi_amp = myamp
                myamp = np.mean((lo_amp, myamp))
            elif Vm_vec_py.max() < targetV:
                lo_amp = myamp
                myamp = np.mean((myamp, hi_amp))

            (axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
                run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, \
                           axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
            Vm_vec_py = np.array(Vm_vec)

        print "D=%.2fum, r=%.2fmm: %.2fnA or %.2fmA extra-cellular stimulation, max change in Vm=%.2fmV" % \
                (fiber_diam, elec_dist, dummy_stim.amp, dummy_stim.amp/10**6, Vm_vec_py.max()-h.v_init)
        threshCur.append(dummy_stim.amp)

    # plot threshold current
    plt.figure()
    plt.plot(D, np.abs(threshCur)/10**6, '*-')
    plt.xlabel("D (um)")
    plt.ylabel("Threshold current magnitude (mA)")
    plt.ylim((0, 1))
    plt.title("D vs. |I|")
    plt.show(block=False)
    plt.savefig('diameter.png', bbox_inches='tight')

    return (axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
# run p1_2
#(axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
#        p1_2(axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
print ""

def p1_3(axons, dummy_elec, dummy_stim, Vm_vec, t_vec):
    """
    Fix fiber diameter at D=11um, use electrode to axon distance 
    ranging from [100um, 200um, 500um, 1mm, 2mm, 5mm, 1cm, 2cm, 5cm], which is
    [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]mm

    Increase cathodic current until Vm=20mV. Record and plot the D vs. current
    relationship
    """
    threshCur = []
    r = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]    # elec_dist in mm
    for rr in r:
        fiber_diam = 11.0  # starting d, in [um]
        num_nodes = 51     # must be int
        my_nseg = 1        # must be int
        elec_dist = rr     # starting r, in [mm]
                
        lo_amp = -10**2       # low cathodic current limit [nA]
        hi_amp = -8*10**10    # high cathodic current  limit [nA]
        
        # starting in the middle -- binary search
        myamp = np.mean((lo_amp, hi_amp))      
        
        (axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
                run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, \
                           axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
        Vm_vec_py = np.array(Vm_vec)
        targetV = h.v_init+20
        
        while np.abs(Vm_vec_py.max() - targetV) > 0.1:
            if Vm_vec_py.max() > targetV:
                hi_amp = myamp
                myamp = np.mean((lo_amp, myamp))
            elif Vm_vec_py.max() < targetV:
                lo_amp = myamp
                myamp = np.mean((myamp, hi_amp))

            (axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
                run_Neuron(fiber_diam, num_nodes, my_nseg, elec_dist, myamp, \
                           axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
            Vm_vec_py = np.array(Vm_vec)

        print "D=%.2fum, r=%.2fmm: %.2fnA or %.2fmA extra-cellular stimulation, max change in Vm=%.2fmV" % \
            (fiber_diam, elec_dist, dummy_stim.amp, dummy_stim.amp/10**6, Vm_vec_py.max()-h.v_init)
        threshCur.append(dummy_stim.amp)

    # plot threshold current
    plt.ylabel("Threshold current magnitude (mA)")
    #plt.ylim((0, 10))
    plt.title("r vs. |I|")
    plt.show(block=False)
    plt.savefig('elec_dist.png', bbox_inches='tight')

    plt.figure()
    plt.semilogy(r, np.abs(threshCur)/10**6, '*-')
    plt.xlabel("Electrode to fiber distance (mm)")
    plt.ylabel("Threshold current magnitude (mA)")
    #plt.ylim((0, 10))
    plt.title("r vs. log(|I|)")
    plt.show(block=False)
    plt.savefig('elec_dist_log.png', bbox_inches='tight')

    return (axons, dummy_elec, dummy_stim, Vm_vec, t_vec)
# run p1_3
(axons, dummy_elec, dummy_stim, Vm_vec, t_vec) = \
        p1_3(axons, dummy_elec, dummy_stim, Vm_vec, t_vec)

def p1_1():
    Ra = []
    Cm = []
    Rm = []
    D = []
    for i in np.arange(1,21+1,1):
        params = calcParams(i, 51, 1, 1, 1)
        D.append(i)
        Ra.append(params.my_Ra) # [ohm]
        Cm.append(params.my_Cm) # [uF]
        Rm.append(params.my_Rm) # [ohm]

    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].plot(D, np.array(Ra)/(10**6), '-*')  # convert Ohm to MOhm
    axarr[0].set_title('Internodal Axoplasmic resistance')
    axarr[0].set_ylabel('MOhm')

    axarr[1].plot(D, np.array(Cm)/(10**-6), '-*') # convert uF to pF
    axarr[1].set_title('Nodal Membrane capacitance')
    axarr[1].set_ylabel('pF')

    axarr[2].plot(D, np.array(Rm)/(10**6), '-*') # convert Ohm to MOhm
    axarr[2].set_title('Nodal Membrane resistance')
    axarr[2].set_ylabel('MOhm')
    axarr[2].set_xlabel('D (um)')

    plt.show(block=False)
    plt.savefig('p1_1.png', bbox_inches='tight')
    return (Ra,Cm,Rm,D)
Ra,Cm,Rm,D = p1_1()
