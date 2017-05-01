from PNS import *

def calcTotalNodes():
    """
    For a fiber, there must be at least 20 nodes before
    the cathodic electrode, and 20 nodes after the cathodic
    electrode.

    Max distance between the two electrodes is 30mm
    """
    testBunch = Bunch()
    testBunch.fiber_diam = 1.0 # [um]
    testBunch.scale_diam = 0.7
    testBunch.num_nodes = 51
    testBunch.num_seg = 1
    testBunch.rho_a = 100.0     # [ohm-cm] 
    testBunch.rho_e = 500       # [ohm-cm]
    testBunch.node_cm = 1.0     # [uF/cm^2]
    testBunch.node_len = 1.5    # [um]
    testBunch.scale_len = 100.0

    calc = MyelinatedNerve.calcParams(testBunch)
    internodal_len = calc.internodal_len
    total_nodes = 20+20+np.ceil(30.0*1000/internodal_len)
    print "total_nodes=", total_nodes
    return int(total_nodes)

# Initialize sweeney nerve
initBunch = Bunch()
initBunch.fiber_diam = 10.0 # [um]
initBunch.scale_diam = 0.7
initBunch.num_nodes = calcTotalNodes()
initBunch.num_seg = 1
initBunch.rho_a = 100.0      # [ohm-cm]
initBunch.rho_e = 500.0       # [ohm-cm]
initBunch.node_cm = 1.0     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0
nerve = MyelinatedNerve(initBunch, sweeneyCh)  

# Initialize stimulating bipolar electrodes
myamp = -0.15e-3/1e-9
mydur = 2 # ms
mydel = 1 # delay of stimulus, in ms
my_dsep = 1    # [mm] separation between electrodes
electrode = DummyBipolarElectrode(mydel, myamp, mydur, my_dsep)
#electrode = DummyElectrode(mydel, myamp, mydur)

# Initialize simulation
simBunch = Bunch()
simBunch.dt = 0.005     # [ms]
simBunch.tstop = 10.0   # [ms]
simBunch.v_init = -80   # [mV]
simBunch.celsius = 37   # [deg C]
simBunch.elec_node_idx = 20 
simBunch.elec_dist = 1  # [mm]
simBunch.node_offset = 0 # [um] offset1 in calcDist
sim = StimSim(nerve, electrode, simBunch)   # active Sweeney

"""
sim.run_sim()
Vm_vec_py = np.array(sim.Vm_vec)
t_vec_py = np.array(sim.t_vec)
plt.figure()
for i in np.arange(1.0, 1.3, 0.005):
    #v_alpha = (i-1.0)/(0.1)
    plt.plot(np.arange(nerve.params.num_nodes), Vm_vec_py[:, int(i/h.dt)], 'r-')
plt.xlabel('node')
plt.ylabel('mV')
plt.ylim((-90, 20))
plt.show(block=False)
"""

def bipolar_effects(fiber_diam):
    durations = [0.02]
    elec_dist = [0.25, 0.5, 0.75, 1.0]
    node_offset = [0, 0.25, 0.5]
    dsep = [1, 5, 10, 15, 20]

    lo_lim = -0.01e-3/1e-9        # cathodic current low limit [nA]
    hi_lim = -10e-3/1e-9         # cathodic current high limit [nA]

    threshCur = np.empty([len(node_offset), len(elec_dist)])
    threshCur_low = np.empty([len(node_offset), len(elec_dist)])
    threshCur_hi = np.empty([len(node_offset), len(elec_dist)])
    offset_cur = np.empty([len(dsep), len(elec_dist)])
    
    sim.nerve.set_fiber_diam(fiber_diam)

    for i in np.arange(len(node_offset)):
        simBunch.node_offset = node_offset[i]*(nerve.params.internodal_len)
        offset_cur[:] = np.nan
        for j in np.arange(len(elec_dist)):
            sim.simParams.elec_dist = elec_dist[j]
            for k in np.arange(len(dsep)):
                electrode.dsep = dsep[k]
                s = strength_duration(sim, durations, lo_lim, hi_lim, AP_status, AP_msgFn)
                offset_cur[k, j] = np.abs(s[:]) # should be only one value here
        threshCur[i,:] = np.mean(offset_cur, 0)
        threshCur_low[i,:] = np.min(offset_cur, 0)
        threshCur_hi[i,:] = np.max(offset_cur, 0)

    threshCur = threshCur/1e6
    threshCur_low = threshCur_low/1e6
    threshCur_hi = threshCur_hi/1e6

    plt.figure()
    for i in np.arange(len(node_offset)):
        yerr = [threshCur_low[i,:]-threshCur[i,:], threshCur_hi[i,:]-threshCur[i,:]]
        plt.errorbar(elec_dist, threshCur[i,:], yerr=yerr, fmt='--o', label='node_offset=%.2f' %(node_offset[i]))
    plt.xlabel('elec_dist [mm]')
    plt.ylabel('threshCur (mA)')
    plt.legend()
    plt.title("%d um: bipolar effects" % (fiber_diam))
    plt.show(block=False)

    return threshCur, threshCur_low, threshCur_hi

def test_one_fiber(fiber_diam):
    durations = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.2, 1.5]  # PW [ms]
    elec_dist = [0.5, 0.75, 1.0, 1.25, 1.5]                     # vertical dist [mm]
    node_offset = [0, 0.25, 0.5]                                # frac of internodal len
    dsep = [1, 5, 10, 15, 20]                                   # lateral elec sep [mm]

    lo_lim = -0.05e-3/1e-9  # cathodic current low limit [nA]
    hi_lim = -10e-3/1e-9    # cathodic current high limit [nA]

    sim.nerve.set_fiber_diam(fiber_diam)
    threshCur = np.empty([len(elec_dist), len(node_offset), len(dsep), len(durations)])
    
    for i in np.arange(len(elec_dist)):
        sim.simParams.elec_dist = elec_dist[i]
        for j in np.arange(len(node_offset)):
            simBunch.node_offset = node_offset[j]*(nerve.params.internodal_len)
            for k in np.arange(len(dsep)):
                electrode.dsep = dsep[k]
                s = strength_duration(sim, durations, lo_lim, hi_lim, AP_status, AP_msgFn)
                threshCur[i, j, k] = s[:]

    print "Saving result for %d um" % (fiber_diam)
    np.save('threshCur_bipolarStim_%dum.npy' % (fiber_diam), threshCur)
    return threshCur

def monopolar_effects(fiber_diam):
# use one electrode
    durations = [0.02]
    elec_dist = [1, 1.5, 2]
    node_offset = [0, 0.25, 0.5]

    lo_lim = -0.05e-3/1e-9        # cathodic current low limit [nA]
    hi_lim = -10e-3/1e-9         # cathodic current high limit [nA]

    threshCur = []
    sim.nerve.set_fiber_diam(fiber_diam)
    for dist in elec_dist:
        sim.simParams.elec_dist = dist  # change elec vertical distance
        for offset in node_offset:
            simBunch.node_offset = offset*(nerve.params.internodal_len) # change offset
            threshCur.extend(strength_duration(sim, durations, lo_lim, hi_lim, AP_status, AP_msgFn))

    threshCur = np.array(threshCur)/1e6
    k = 0
    plt.figure()
    for i in np.arange(len(node_offset)):
        plt.plot(elec_dist, threshCur[k:k+len(elec_dist)], '-*', label='node_offset=%.2f' % (node_offset[i]))
        k = k+len(node_offset)
    plt.xlabel('elec_dist [mm]')
    plt.ylabel('range of threshCur (mA)')
    plt.legend()
    plt.title('%d um, monopolar effects' % (fiber_diam))
    plt.show(block=False)

    return threshCur

"""
# --- Try to get an action potential
durations = [0.1] #[0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.2, 1.5] # [ms]
node_offset = [0, 0.25, 0.5]    # fraction of internodal len [um] 
#dsep = [1, 5, 10, 15, 20]

lo_lim = -0.05e-3/1e-9        # cathodic current low limit [nA]
hi_lim = -5e-3/1e-9         # cathodic current high limit [nA]

threshCur = []
for offset in node_offset:
    #for sep in dsep:
    #    electrode.dsep = sep
        simBunch.node_offset = offset*(nerve.params.internodal_len) # convert to [um]
        threshCur.append(strength_duration(sim, durations, lo_lim, hi_lim, AP_status, AP_msgFn))
"""

"""
threshCur = (strength_duration(sim, durations, lo_lim, hi_lim, AP_status, AP_msgFn))

plt.figure()
plt.plot(durations, np.abs(np.array(threshCur))/1e6, '-*')
plt.xlabel('durations (ms)')
plt.ylabel('Current amplitude (mA)')
plt.title('strength-duration relationship for %dum fibers' % nerve.params.fiber_diam)
#plt.legend([str(i)+'ms' for i in durations])
plt.show(block=False)
"""
