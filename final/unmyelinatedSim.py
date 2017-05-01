from PNS import *


#--- Create unmyelianted axon first
# Unmyelinated axon initialization parameter bunch
initBunch = Bunch()
initBunch.fiber_diam = 0.5    #[um]
initBunch.rho_a = 100.0     # axioplasmic resistivity [ohm-cm]
initBunch.rho_e = 500.0     # extracellular resistivity [ohm-cm]
initBunch.cm = 1.0          # [uF/cm^2] membrane cap
initBunch.L = 1000          # [um] ~the length of 10 nodes of 10um diameter mylinated axon
initBunch.segLen = 2.5      # [um] limits length of each segment
initBunch.nseg = makeOdd((initBunch.L/initBunch.segLen))   # always odd

smallFiber = UnmyelinatedNerve(initBunch, sweeneyCh)

#--- Stimulation electrode
myAmp = -5.25e-3/1e-9    # amplitude [nA]
myDel = 5.0             # start at t=5ms [ms]
myDur = 1             # duration [ms]
electrode = DummyElectrode(myDel, myAmp, myDur)

#--- Unmyelinated simulation initialization parameter bunch for testing
simBunch = Bunch()
simBunch.dt = 0.005       # [ms]
simBunch.v_init = -70.0   # [mV]
simBunch.celsius = 37.0   # [deg C]
simBunch.elec_dist = 0.25  # [mm]
simBunch.elecLocX = 0.5    # 0 to 1, normalized location along unmeylinated axon
simBunch.tstop = 10.0     # [ms]

sim = StimSim(smallFiber, electrode, simBunch)

#--- simulation
#sim.run_sim()
#Vm_vec_py = np.array(sim.Vm_vec)
#t_vec_py = np.array(sim.t_vec)

#--- Plotting utilities
def plot_segV(Vm_vec_py, t_vec_py, segIdx):
# Plot the time course of voltage at a segment
    if segIdx >= Vm_vec_py.shape[0]:
        print "Only %d segments, %d invalid!" % (Vm_vec_py.shape[0], segIdx)
        return
    plt.figure()
    plt.plot(t_vec, Vm_vec_py[segIdx, :])
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.title('Time course @ %.2fum' % (smallFiber.idx2len(segIdx)))
    plt.show(block=False)

def plot_segV_together(Vm_vec_py, t_vec_py, listIdx):
    plt.figure()
    for i in range(len(listIdx)):
        v_alpha = (i-0.0)/(len(listIdx)) + 0.1
        plt.plot(t_vec_py, Vm_vec_py[listIdx[i],:], 'r-', alpha=v_alpha)
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.title('Time course for %.2fum to %.2fum' % \
             (smallFiber.idx2len(listIdx[0]), smallFiber.idx2len(listIdx[-1])))
    plt.show(block=False)

def plot_allV(Vm_vec_py, tStart=4.0, tEnd=10.0, dt=0.1):
# Plot voltage over entire axon, different
# curves for different time periods
    if tStart > h.tstop:
        print "tStart > tstop!"
        return
    if tStart >= tEnd:
        print "tStart must be less than tEnd!"
        return

    plt.figure()
    L = np.array([seg.x for seg in smallFiber.axon]) * smallFiber.axon.L
    for i in np.arange(tStart, tEnd, dt):
        # later is more opaque
        #import pdb; pdb.set_trace()
        v_alpha = (i-tStart)/(tEnd-tStart)
        #plt.plot(L, Vm_vec_py[:, int(i/h.dt)], 'r-', alpha=v_alpha)
        plt.plot(L, Vm_vec_py[:, int(i/h.dt)])
    plt.xlabel('length (um)')
    plt.ylabel('mV')
    plt.show(block=False)

#plot_allV(Vm_vec_py, tStart=4.5, tEnd=8.5, dt=0.05)
#plot_segV_together(Vm_vec_py, t_vec_py, [smallFiber.len2Idx(l) for l in np.arange(100,1000,100)])
#plt.xlim((4,4+electrode.dummy_stim.dur+1))
#plt.show(block=False)

# --- Actual work for report
electrode.setDelay(2)
durations = [1, 2.5, 3, 5, 10]
diameters = [0.2, 0.5, 0.7, 1.0, 1.2, 1.5]
lo_lim = -0.5e-3/1e-9        # cathodic current low limit [nA]
hi_lim = -10e-3/1e-9         # cathodic current high limit [nA]
threshCur = []
sim.simParams.elec_dist = 0.5  # 0.25mm
print "elec_dist=%.2f mm" % (sim.simParams.elec_dist)

for i in range(len(durations)):
    print "PW=%.2fms: " % (durations[i])
    electrode.setDur(durations[i])
    sim.change_tstop(np.max([10, 5+durations[i]]))
    threshCur.append(strength_diam(sim, diameters, lo_lim, hi_lim, AP_status, AP_msgFn))


