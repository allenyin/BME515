from PNS import *

# Initialize sweeney nerve
initBunch = Bunch()
initBunch.fiber_diam = 12.0 # [um]
initBunch.scale_diam = 0.7
initBunch.num_nodes = 51
initBunch.num_seg = 1
initBunch.rho_a = 100.0      # [ohm-cm]
initBunch.rho_e = 500.0       # [ohm-cm]
initBunch.node_cm = 1.0     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0
nerve = MyelinatedNerve(initBunch, sweeneyCh)  

# Initialize stimulating electrode
myamp = -0.115e-3/1e-9
mydur = 0.2 # ms
mydel = 5
electrode = DummyElectrode(mydel, myamp, mydur)

# Initialize simulation
simBunch = Bunch()
simBunch.dt = 0.005     # [ms]
simBunch.tstop = 10.0   # [ms]
simBunch.v_init = -80   # [mV]
simBunch.celsius = 37   # [deg C]
simBunch.elec_node_idx = initBunch.num_nodes/2
simBunch.elec_dist = 1  # [mm]
sim = StimSim(nerve, electrode, simBunch)   # active Sweeney

# --- Actual work for report
electrode.setDelay(2)
durations = [0.2, 0.5, 0.7, 1.0, 2.5, 3.0, 5.0, 10.0]
diameters = [1, 1.2, 1.5, 2.2, 2.5, 3]
lo_lim = -0.1e-3/1e-9        # cathodic current low limit [nA]
hi_lim = -10e-3/1e-9         # cathodic current high limit [nA]
threshCur = []
sim.simParams.elec_dist = 1.386  # [mm] 
print "elec_dist=%.2f mm" % (sim.simParams.elec_dist)

for i in range(len(durations)):
    print "PW=%.2fms: " % (durations[i])
    electrode.setDur(durations[i])
    sim.change_tstop(np.max([10, 5+durations[i]]))
    threshCur.append(strength_diam(sim, diameters, lo_lim, hi_lim, AP_status, AP_msgFn))

plt.figure()
for i in range(len(durations)):
    plt.plot(diameters, np.abs(np.array(threshCur[i]))/1e6)
plt.xlabel('diameters (um)')
plt.ylabel('Current amplitude (mA)')
plt.title('diameter-current relationship for B-fibers')
plt.legend([str(i)+'ms' for i in durations])
plt.show(block=False)

plt.axhline(y=3.5, color='r', linewidth=4)
plt.show(block=False)
plt.savefig('B_fibers_diam_current.png', bbox_inches='tight')
