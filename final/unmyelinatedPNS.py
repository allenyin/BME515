import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from neuron import h, gui
from bunch import Bunch
np.seterr(all='raise')

def makeOdd(n):
    if n%2==0:
        return int(n+1)

def initParams(fiber_diam, rho_a=100.0, rho_e=500.0, cm=1.0, L=1000, segLen=2.5):
    """
    Return initBunch based on fiber_diam
    """
    initBunch = Bunch()
    initBunch.fiber_diam = fiber_diam # [um]
    initBunch.rho_a = rho_a   # axioplasmic resistivity [ohm-cm]
    initBunch.rho_e = rho_e   # extracellular resistivity [ohm-cm]
    initBunch.cm = cm         # [uF/cm^2] membrane cap
    initBunch.L = L           # [um] ~the length of 10 nodes of 10um diameter mylinated axon
    initBunch.segLen = segLen # [um] limits length of each segment
    initBunch.nseg = makeOdd((L/segLen))   # always odd
    return initBunch

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

def createElec(delay, amp, dur):
    """
    Return a Neuron electrode section with the given properties,
    as well as a stimulation IClamp
    """
    dummy_elec = h.Section('dummy_elec')
    dummy_stim = h.IClamp(dummy_elec(0.5))
    dummy_stim.delay = delay    # start after this time [ms]
    dummy_stim.amp = amp        # amplitude [nA]. Cathodic current is negative
    dummy_stim.dur = dur        # duration, aka pulsewidth [ms]
    return dummy_elec, dummy_stim

def calcDist(axon, elecLoc, elec_dist_um):
    """
    elecLoc is the location of the electrode in terms of
    normalized axon location (0 to 1).

    Returns a vector of distances of each of the axon's
    segment to the electrode
    """
    distToElec = np.array([seg.x-elecLoc for seg in axon]) * axon.L
    return np.sqrt(np.abs(distToElec)**2 + elec_dist_um**2)


