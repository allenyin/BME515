import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from neuron import h, gui
from bunch import Bunch
np.seterr(all='raise')

'''
Test if sweeney on action potential generation for 
unmyelinated axons.
'''
def sweenyCh(axon):
    axon.insert('sweeney')
    axon.ena = 35.64    # [mV]
    return axon

def createAxon(initBunch, chFunc=sweenyCh):
    axon = h.Section('axon')
    axon.push()
    axon.nseg = initBunch.nseg
    axon.diam = initBunch.fiber_diam
    axon.L = initBunch.L
    axon.Ra = initBunch.rho_a
    axon.cm = initBunch.cm
    axon.insert('extracellular')
    axon = chFunc(axon)
    return axon

def calcDist(axon, elecLoc, elec_dist_um):
    """
    elecLoc is the location of the electrode in terms of
    normalized axon location (0 to 1).

    Returns a vector of distances of each of the axon's
    segment to the electrode
    """
    distToElec = np.array([seg.x-elecLoc for seg in axon]) * axon.L
    return np.sqrt(np.abs(distToElec)**2 + elec_dist_um**2)

def makeOdd(n):
    if n%2==0:
        return int(n+1)

# Let's use group-C fibers, so diameter is 0.2-1.5um
initBunch = Bunch()
initBunch.fiber_diam = 1.5  # [um]
initBunch.rho_a = 100.0   # axioplasmic resistivity [ohm-cm]
initBunch.rho_e = 500.0   # extracellular resistivity [ohm-cm]
initBunch.cm = 1.0   # [uF/cm^2] membrane cap
initBunch.L = 1000       # [um] ~the length of 10 nodes of 10um diameter mylinated axon
initBunch.segLen = 2.5    # [um] limits length of each segment
initBunch.nseg = makeOdd((initBunch.L/initBunch.segLen))   # always odd

# instrumentation params
dummy_elec = h.Section('dummy_elec')
dummy_stim = h.IClamp(dummy_elec(0.5))
dummy_stim.delay = 5.0  # start at t=5ms [ms]
dummy_stim.amp = -25.0e-3/1e-9     # amplitude [nA]. Cathodic current is negative
dummy_stim.dur = 2.5    # duration, aka pulsewidth [ms] -- 0.1ms=100us

# simulation parameters
h.dt = 0.005      # [ms]
h.tstop = np.max([10,5+dummy_stim.dur])      # [ms]
h.v_init = -80.0  # [mV]..sweeney has resting potential=-80mV
h.celsius = 37.0  # [deg C]
num_timesteps = int(h.tstop/h.dt)

# electrode params
elec_dist = 0.6                 # shortest distance from electrode to fiber [mm]
elec_dist_um = elec_dist*10**3  # distance in [um]

# ---- Simulation parts
smallFiber = createAxon(initBunch)
rVec = calcDist(smallFiber, 0.2, elec_dist_um)  # elec at 0.2 from the start of fiber
segs = list(smallFiber.allseg())

# record Vm(t) at all nodes
Vm_vec = []
for seg in smallFiber:
    Vm_vec.append(h.Vector(num_timesteps, 0))
    Vm_vec[-1].record(seg._ref_v)

# time vector
t_vec = h.Vector()
t_vec.record(h._ref_t)

#------- Simulation control
def run_sim():
    h.finitialize(h.v_init)
    i = 0
    while (h.t < h.tstop):
        for seg in smallFiber:
            #print i, ",", seg.x
            #if np.abs(dummy_stim.i)>0:
            seg.e_extracellular = (10**-2)*initBunch.rho_e*dummy_stim.i/(4*pi*rVec[i])
            i = i+1
        i = 0
        h.fadvance()
    Vm_vec_py = np.array(Vm_vec)
    return Vm_vec_py

Vm_vec_py = run_sim()

#--- Plotting
def len2Idx(x):
# Convert from length to index of the node at that location
    L = initBunch.L
    segPos = np.array([seg.x for seg in smallFiber])*L
    return np.argmin(np.abs(segPos-x))
    
def frac2Idx(frac):
# Convert fraction of length (0 to 1) to index of the node
# at that location
    segPos = np.array([seg.x for seg in smallFiber])
    return np.argmin(np.abs(segPos-frac))
    
def idx2len(idx):
    L = initBunch.L
    segPos = np.array([seg.x for seg in smallFiber])*L
    return segPos[idx]

def plot_segV(Vm_vec_py, segIdx):
# Plot the time course of voltage at a segment
    if segIdx >= initBunch.nseg:
        print "Only %d segments, %d invalid!" % (initBunch.nseg, segIdx)
        return
    plt.figure()
    plt.plot(t_vec, Vm_vec_py[segIdx, :])
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.title('Time course @ %.2fum' % (idx2len(segIdx)))
    plt.show(block=False)

def plot_segV_together(Vm_vec_py, listIdx):
    plt.figure()
    for i in range(len(listIdx)):
        v_alpha = (i-0.0)/(len(listIdx)) + 0.1
        plt.plot(t_vec, Vm_vec_py[listIdx[i],:], 'r-', alpha=v_alpha)
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.title('Time course for %.2fum to %.2fum'%(idx2len(listIdx[0]), idx2len(listIdx[-1])))
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
    L = np.array([seg.x for seg in smallFiber]) * initBunch.L
    for i in np.arange(tStart, tEnd, dt):
        # later is more opaque
        #import pdb; pdb.set_trace()
        v_alpha = (i-tStart)/(tEnd-tStart)
        #plt.plot(L, Vm_vec_py[:, int(i/h.dt)], 'r-', alpha=v_alpha)
        plt.plot(L, Vm_vec_py[:, int(i/h.dt)])
    plt.xlabel('length (um)')
    plt.ylabel('mV')
    plt.show(block=False)

'''
plot_segV(Vm_vec_py, 50)
plt.savefig('1.png', bbox_inches='tight')
plot_segV(Vm_vec_py, 80)
plt.savefig('2.png', bbox_inches='tight')
plot_segV(Vm_vec_py, 200)
plt.savefig('3.png', bbox_inches='tight')
'''

plot_allV(Vm_vec_py, tStart=4.5, tEnd=8.5, dt=0.05)
#plt.savefig('all.png', bbox_inches='tight')

plot_segV_together(Vm_vec_py, [len2Idx(l) for l in np.arange(100,1000,100)])
plt.xlim((4,4+dummy_stim.dur+1))
plt.show(block=False)
#plt.savefig('timeCourseSegs.png', bbox_inches='tight')
