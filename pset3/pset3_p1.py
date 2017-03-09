from PNS import *
import numpy as np
import pdb


# Initialize sweeney nerve
initBunch = Bunch()
initBunch.fiber_diam = 12.0 # [um]
initBunch.scale_diam = 0.6
initBunch.num_nodes = 51
initBunch.num_seg = 1
initBunch.rho_a = 54.7      # [ohm-cm]
initBunch.rho_e = 500       # [ohm-cm]
initBunch.node_cm = 2.5     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0
nerve = MyelinatedNerve(initBunch, sweeneyCh)  

# Initialize stimulating electrode
myamp = -0.15e-3/1e-9   # 0.15mA cathodic current for 0.2ms duration for 12um fiber
mydur = 0.4 # ms
mydel = 5
elec = DummyElectrode(mydel, myamp, mydur)

# Initialize simulation
simBunch = Bunch()
simBunch.dt = 0.005     # [ms]
simBunch.tstop = 10.0   # [ms]
simBunch.v_init = -80   # [mV]
simBunch.celsius = 37   # [deg C]
simBunch.elec_node_idx = initBunch.num_nodes/2
simBunch.elec_dist = 1  # [mm]
sim1 = StimSim(nerve, elec, simBunch)   # active Sweeney

# Run simulation for different fiber diameter
D = np.arange(2,22,2)   # in um
ampList = np.ones(D.shape)
#ampList = ampList * -1.5e-3/1e-9
ampList[0:3] = [-1e-3/1e-9, -0.5e-3/1e-9, -0.5e-3/1e-9] # roughly current needed to start AP
ampList[3:] = -0.2e-3/1e-9
Vm_list = []
Im_list = []

for (d,a) in zip(D, ampList):
    #pdb.set_trace()
    nerve.set_fiber_diam(d)
    #print "After changing nerve, diam=", nerve.axons[0].diam
    sim1.elec.setAmp(a)
    #print "In sim, diam=", sim1.nerve.axons[0].diam
    sim1.run_sim()
    Vm_list.append(np.array(sim1.Vm_vec))
    Im_list.append(sim1.getImVec())

t_vec_py = np.array(sim1.t_vec)

'''
Plot 1: plot action potential Vm, Im at node 25 and node 45 for all time.
        Blue=node 26, Red=node 40, Green=node 50
        Increasing transparency=Increasing diameter
'''
f, (ax1, ax2) = plt.subplots(2, sharex=True)
for i in range(len(D)):
    alpha = 0.3+i*(0.7/len(D))
    alpha = 1
    # Vm first
    p1 = ax1.plot(t_vec_py[1000:1250], Vm_list[i][26-1, 1000:1250], 'b-', alpha=alpha)
    p2 = ax1.plot(t_vec_py[1000:1250], Vm_list[i][40-1, 1000:1250], 'r-', alpha=alpha)
    p3 = ax1.plot(t_vec_py[1000:1250], Vm_list[i][50-1, 1000:1250], 'g-', alpha=alpha)
    # Im
    ax2.plot(t_vec_py[1000:1250], Im_list[i][26-1, 1000:1250], 'b-', alpha=alpha)
    ax2.plot(t_vec_py[1000:1250], Im_list[i][40-1, 1000:1250], 'r-', alpha=alpha)
    ax2.plot(t_vec_py[1000:1250], Im_list[i][50-1, 1000:1250], 'g-', alpha=alpha)

ax1.set_title('Cathodic stimulation Vm')
ax2.set_title('Cathodic stimulation Im')
ax1.set_ylabel('Vm (mV)')
ax2.set_ylabel('Im (nA)')
ax2.set_xlabel('Time (mS)')

ax1.legend([p1[0],p2[0],p3[0]], ['Center node (node 26)', 'Proximal node (node 40)', 'Distal node (node 50)'])
plt.show(block=False)
f.set_size_inches([11,8])
plt.savefig('p1_node_Vm_Im.png', bbox_inches='tight')

'''
Plot 2, 3: 
Stimulation starts at 5ms, lasts 0.4ms, dt=0.005ms.
So plot Vm, Im @ [5ms, 5.05ms, 5.10ms, 5.15ms, 5.20ms] for all nodes
plot action potential Vm, Im at all nodes at time.

The time is [1000, 1100, 1200, 1300, 1400]
'''
d_idx = [1,3,5,7,9]
f, axarr = plt.subplots(len(d_idx), sharex=True)
for i in np.arange(len(axarr)):
    axarr[i].set_title('Diameter=%d over entire fiber' % (D[d_idx[i]]))
    axarr[i].set_ylabel('Vm (mV)')
    for t in [1000, 1050, 1100, 1150, 1200]:
        axarr[i].plot(np.arange(initBunch.num_nodes), Vm_list[i][:, t])
axarr[-1].set_xlabel('Node')
plt.show(block=False)
plt.savefig('p1_allNode_Vm.png', bbox_inches='tight')
f.set_size_inches([11,11])

f, axarr = plt.subplots(len(d_idx), sharex=True)
for i in np.arange(len(axarr)):
    axarr[i].set_title('Diameter=%d over entire fiber' % (D[d_idx[i]]))
    axarr[i].set_ylabel('Im (nA)')
    for t in [1000, 1050, 1100, 1150, 1200]:
        axarr[i].plot(np.arange(initBunch.num_nodes), Im_list[i][:, t])
axarr[-1].set_xlabel('Node')
plt.show(block=False)
plt.savefig('p1_allNode_Im.png', bbox_inches='tight')
f.set_size_inches([11,11])

'''
Plot 3: plot the peak Im as a function of diameter.
        Take the peak current as the max of the current observed in 
        nodes 40 to 49, as the Im there changes after Im change in the 
        location of stimulation.

        positive Im is measured as current going out
'''
plt.figure()
maxCur = [np.max(I[39:50]) for I in Im_list]
l1, = plt.semilogx(D, maxCur, 'b*-')
plt.xlabel('Fiber diameter (um)')
plt.ylabel('Max transmembrane current (nA)')
plt.title('Fiber diameter vs. Peak Action Potential Current')
plt.show(block=False)
plt.savefig('p1_peak_Im.png', bbox_inches='tight')

