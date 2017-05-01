import numpy as np
pi = np.pi
from matplotlib import pyplot as plt
from neuron import h, gui
from bunch import Bunch
np.seterr(all='raise')

"""
Include the classes and helper methods to create
a myelinated nerve, external mono-polar stimulation electrode,
and run simulation on it.
"""
def makeOdd(n):
    if n%2==0:
        return int(n+1)

"""
# Example Myelinated axon initialization parameter bunch for testing
initBunch = Bunch()
initBunch.fiber_diam = 12.0 # [um]
initBunch.scale_diam = 0.7
initBunch.num_nodes = 51
initBunch.num_seg = 1
initBunch.rho_a = 100.0     # [ohm-cm] 
initBunch.rho_e = 500       # [ohm-cm]
initBunch.node_cm = 1.0     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0

# Myelinated simulation initialization parameter bunch for testing
simBunch = Bunch()
simBunch.dt = 0.005     # [ms]
simBunch.tstop = 10.0   # [ms]
simBunch.v_init = -70.0 # [mV]
simBunch.celsius = 37.0 # [deg C]
simBunch.elec_node_idx = initBunch.num_nodes/2  # for myelinated axon
simBunch.elec_dist = 1.0  # [mm]
"""

# Example unmyelinated axon initialization parameter bunch
initBunch = Bunch()
initBunch.fiber_diam = 1.5  #[um]
initBunch.rho_a = 100.0     # axioplasmic resistivity [ohm-cm]
initBunch.rho_e = 500.0     # extracellular resistivity [ohm-cm]
initBunch.cm = 1.0          # [uF/cm^2] membrane cap
initBunch.L = 1000          # [um] ~the length of 10 nodes of 10um diameter mylinated axon
initBunch.segLen = 2.5      # [um] limits length of each segment
initBunch.nseg = makeOdd((initBunch.L/initBunch.segLen))   # always odd

# Unmyelinated simulation initialization parameter bunch for testing
simBunch = Bunch()
simBunch.dt = 0.005     # [ms]
simBunch.tstop = 10.0   # [ms]
simBunch.v_init = -70.0 # [mV]
simBunch.celsius = 37.0 # [deg C]
simBunch.elec_dist = 1.0  # [mm]
simBunch.elecLoc = 0.2  # 0 to 1, normalized location along unmeylinated axon

def passiveCh(axon):
    node_rm = 2000  # membrance specific resistance [ohm-cm**2] 
    axon.insert('pas')
    axon.g_pas = 1.0/node_rm
    return axon

def passiveSweeneyCh(axon):
    axon.insert('pas')
    axon.g_pas = 0.128 # mho/cm^2, Sweeney

    return axon

def sweeneyCh(axon):
    axon.insert('sweeney')
    axon.ena = 35.64 # [mV]
    return axon

class MyelinatedNerve(object):
    """
    Class to describe a myelinated nerve made of connected axons.
    
    Essentially a wrapper class around a list of axons, each axon
    is a Section in Neuron.

    On initialization, requires parameter Bunch object, and a function
    that sets the channel properties of each axon given an axon.

    parameter Bunch contains the entries:
        fiber_diam:  Fiber diameter [um]
        scale_diam:  Scaling used to obtain node_diam from fiber_diam
        num_nodes:   Number of axon sections.
        num_seg:     Number of segments per axon
        rho_a:       Axoplasmic resistivity [ohm-cm]
        rho_e:       Extracellular resistivity [ohm-cm]
        node_cm:     Nodal membrance capacitance [uF/cm^2]
        node_len:    Node length [um]
        scale_len:   Scaling used to obtain myelin length from fiber_diam

    Each instance contains the following variables:
        params: Bunch object with params relevant to the nerve, include given
        given:  Bunch object with params given by user
        chFunc: function used to add channels to an axon section
        axons:  list of hocObjects referring to axon sections
    """
    def __init__(self, paramBunch, chFunc):
        self.initParams(paramBunch)
        self.chFunc = chFunc
        self.initAxons()

    def set_fiber_diam(self, fiber_diam):
        self.given.fiber_diam = fiber_diam
        self.params = self.given.copy()
        self.params.update(self.calcParams(self.given))
        self.params = Bunch(self.params)
        for axon in self.axons:
            axon = self.editOneAxon(axon)

    def add_nodes(self, n_new):
        self.given.num_nodes = self.given.num_nodes + n_new
        self.params = self.given.copy()
        self.params.update(self.calcParams(self.given))
        self.params = Bunch(self.params)
        #import pdb; pdb.set_trace()
        extra_axons = self.createAxons(n_new)
        extra_axons[0].connect(self.axons[-1],1,0)
        self.axons.extend(extra_axons)

    def initAxons(self):
        #import pdb; pdb.set_trace()
        self.axons = self.createAxons(self.params.num_nodes)

    def createAxons(self, n):
        """
        Create n axons
        """
        newAxons = []
        for i in range(n):
            newAxons.append(h.Section('axon'))
            newAxons[i].push()
            newAxons[i] = self.editOneAxon(newAxons[i])
            h.pop_section()
        for i in range(n-2):
            newAxons[i+1].connect(newAxons[i],1,0)
        return newAxons
    
    def editOneAxon(self, axon):
        axon.nseg = self.params.num_seg
        axon.diam = self.params.node_diam
        axon.L = self.params.node_len
        axon.Ra = self.params.rho_a_adj
        axon.cm = self.params.node_cm
        axon.insert('extracellular')
        axon = self.chFunc(axon)
        return axon
       
    def initParams(self, paramBunch):
        self.given = paramBunch
        self.params = paramBunch.copy()
        self.params.update(self.calcParams(self.given))
        self.params = Bunch(self.params)

    def calcDist(self, elec_node_idx, elec_dist_um_y, offset1, offset2):
        """
        Extracellular electrode is placed above the node
        with index = elec_node_idx, at a distance of
        elec_dist_um [um] away.

        Returns vector of distances from each of the axon's
        nodes to the electrode, in [um]
        """
        internodal_len = self.params.internodal_len
        node_dist = np.arange(self.params.num_nodes)*internodal_len # dist of node from head
        e1_dist = elec_node_idx*internodal_len-offset1
        e2_dist = e1_dist + offset2
        lateral_dist1 = np.abs(node_dist-e1_dist)
        lateral_dist2 = np.abs(node_dist-e2_dist)
        return np.sqrt((lateral_dist1**2 + elec_dist_um_y**2)), \
               np.sqrt((lateral_dist2**2 + elec_dist_um_y**2))

    @staticmethod
    def calcParams(given):
        calc = Bunch()

        calc.node_diam = given.scale_diam * given.fiber_diam
        calc.node_diam_cm = calc.node_diam/10**4
        calc.node_len_cm = given.node_len/10**4
        calc.myelin_len = given.scale_len * given.fiber_diam
        calc.myelin_len_cm = calc.myelin_len/10**4

        calc.internodal_len = calc.myelin_len + 0.5*given.node_len
        calc.rho_a_adj = \
           given.rho_a*((calc.node_len_cm + calc.myelin_len_cm)/calc.node_len_cm)
        return calc
   
class UnmyelinatedNerve(object):
    """
    Class to describe an unmyelinated nerve made of connected axons.

    Essentially a wrapper class around a list of axons, each axon is a Section
    in Neuron.

    On initialization, requires parameter Bunch object, and a function that sets
    the channel properties of each axon given an axon.

    parameter Bunch contains entries:
        fiber_diam: Fiber diameter [um]
        rho_a = 100.0     # axioplasmic resistivity [ohm-cm]
        rho_e = 500.0     # extracellular resistivity [ohm-cm]
        cm = 1.0          # [uF/cm^2] membrane cap
        L = 1000          # [um] ~the length of 10 nodes of 10um diameter mylinated axon
        segLen = 2.5      # [um] limits length of each segment
        nseg = makeOdd((initBunch.L/initBunch.segLen))   # always odd

    Each instance contains the following variables:
        params: Bunch object with params relevant to the nerve, include given
        given:  Bunch object with params given by user
        chFunc: function used to add channels to an axon section
        axons:  list of hocObjects referring to axon sections
    """
    def __init__(self, paramBunch, chFunc):
        self.initParams(paramBunch)
        self.chFunc = chFunc
        self.initAxons()

    def set_fiber_diam(self, fiber_diam):
        self.params.fiber_diam = fiber_diam
        self.axon.diam = self.params.fiber_diam

    def initAxons(self):
        self.axon = h.Section('axon')
        self.axon.push()
        self.axon.nseg = self.params.nseg
        self.axon.diam = self.params.fiber_diam
        self.axon.L = self.params.L
        self.axon.Ra = self.params.rho_a
        self.axon.cm = self.params.cm
        self.axon.insert('extracellular')
        self.axon = self.chFunc(self.axon)

    def initParams(self, paramBunch):
        self.params = paramBunch

    def len2Idx(self, l_um):
        """
        Given length [um], convert to index of the segment
        corresponding to that length
        """
        assert l_um < self.params.L
        segPos = np.array([seg.x for seg in self.axon]) * self.params.L
        return np.argmin(np.abs(segPos-l_um))

    def frac2Idx(self, frac):
        """
        Given fraction of the length (0 to 1) along the axon
        to the index of segment corresponding to that location
        """
        assert frac>0 and frac<=1
        segPos = np.array([seg.x for seg in self.axon])
        return np.argmin(np.abs(segPos - frac))

    def idx2len(self, idx):
        """
        Given index of the segment, return the length of the axon there
        """
        assert idx < self.params.nseg
        L = self.params.L
        segPos = np.array([seg.x for seg in self.axon])*L
        return segPos[idx]

    def calcDist(self, elecLocX, elec_dist_um_y):
        """
        elecLocX is the location of the electrode
        along the fiber in normalized location (0 to 1)

        elec_dist_um_y is the shortest perpendicular separation
        between the electrode and the fiber. In [um]

        Return vector of dis from each seg to electorde in [um]
        """
        assert elecLocX > 0 and elecLocX < 1
        distToElec = np.array([seg.x-elecLocX for seg in self.axon]) * self.axon.L
        return np.sqrt(np.abs(distToElec)**2 + elec_dist_um_y**2)

class DummyElectrode(object):
    """
    Class to describe the dummy electrode attached to the axon.
    """

    def __init__(self, delay, amp, dur):
        self.dummy_elec = h.Section('dummy_elec')
        self.dummy_stim = h.IClamp(self.dummy_elec(0.5))
        self.dummy_stim.delay = delay
        self.dummy_stim.amp = amp
        self.dummy_stim.dur = dur

    def setDelay(self, delay):
        self.dummy_stim.delay = delay

    def setAmp(self, amp):
        self.dummy_stim.amp = amp

    def setDur(self, dur):
        self.dummy_stim.dur = dur

class DummyBipolarElectrode(DummyElectrode):
    """
    Class to describe bipolar dummy electrodes attached to the axon,
    with separation distance d_sep

    Essentially, same as DummyElectrode but with extra property d_sep.
    """
    def __init__(self, delay, amp, dur, dsep):
        DummyElectrode.__init__(self, delay, amp, dur)
        self.dsep = dsep

    def getE1_current(self):
        return self.dummy_stim.i

    def getE2_current(self):
        return -self.dummy_stim.i

class StimSim(object):
    """
    Class to describe an extracellular stimulation simulation.
    Requires a MyelinatedNerve or UnmyelinatedNerve and a DummyElectrode object.

    In addition, have options to change simulation environment such as
    time, temperature, resting potential.
    
    simulation initialization parameter bunch example
        simBunch = Bunch()
        simBunch.dt = 0.005     # [ms]
        simBunch.tstop = 10.0   # [ms]
        simBunch.v_init = -70.0 # [mV]
        simBunch.celsius = 37.0 # [deg C]
        simBunch.elec_node_idx = initBunch.num_nodes/2
        simBunch.elec_dist = 1  # [mm]

    For UnmyelinatedNerve, use simBunch.elecLocX instead of simBunch.elec_node_idx.
        simBunch.elecLocX is between 0 and 1, represents normalized distance along
        the nerve fiber.

    Each instance contains:
        nerve: MyelinatedNerve object reference
        elec:  DummyElectrode object reference
        simParams: simulation parameter bunch
        Vm_vec:   Hoc vector containing Vm      (in mV)
        t_vec:    Hoc vector containing time    (in mS)
        Im_vec:   Hoc vector containing transmembrane current. (in mA/cm^2)
                  Initialized if constructor 'recordI' is true
                  
    In this case, the Vm_vec and Im_vec will contain that of all MyelinatedNerves

    Bipolar electrode:
        The field applied at each node is the superposition of the two separated
        point sources. Using superposition of two monopoles is needed since the
        distance from the source to the fiber is not orders of magnitude larger
        than the source separation.

    ** Bug in here: If the nerve increase/decreases length, this bugs out, no
                    time to fix and test this.
    """
    def __init__(self, nerve, elec, simBunch):
        self.nerve = nerve
        self.elec = elec
        self.simParams = simBunch
        num_timesteps = int(simBunch.tstop/simBunch.dt)+1
        self.simParams.num_timesteps = num_timesteps
        self.create_vectors()

    def create_vectors(self):
        self.Vm_vec = []
        self.Im_vec = []
        if type(self.nerve) is MyelinatedNerve:
            for i in range(self.nerve.params.num_nodes):
                self.Vm_vec.append(h.Vector(self.simParams.num_timesteps, 0))
                self.Vm_vec[i].record(self.nerve.axons[i](0.5)._ref_v)
                self.Im_vec.append(h.Vector(self.simParams.num_timesteps, 0))
                self.Im_vec[i].record(self.nerve.axons[i](0.5)._ref_i_membrane)
        elif type(self.nerve) is UnmyelinatedNerve:
            for seg in self.nerve.axon:
                self.Vm_vec.append(h.Vector(self.simParams.num_timesteps, 0))
                self.Vm_vec[-1].record(seg._ref_v)
                self.Im_vec.append(h.Vector(self.simParams.num_timesteps, 0))
                self.Im_vec[-1].record(seg._ref_i_membrane)
        self.t_vec = h.Vector()
        self.t_vec.record(h._ref_t)

    def getImVec(self):
        # convert the Im_vec [mA/cm^2] by multiplying area [um^2], to give nA
        # for single nerve use
        if type(self.nerve) is MyelinatedNerve:
            area = pi * self.nerve.params.node_diam * self.nerve.params.node_len
            print "getImVec: diam=", self.nerve.axons[0].diam
            print "getImVec: area=", area 
            Im_vec_nA = np.array(self.Im_vec)
            Im_vec_nA = Im_vec_nA*area
            return Im_vec_nA
        elif type(self.nerve) is UnmyelinatedNerve:
            print "getImVec() not implemented for UnmyelinatedNerve yet!"
            return

    def change_tstop(self, tstop):
        self.simParams.tstop = tstop
        self.simParams.num_timesteps = int(self.simParams.tstop/self.simParams.dt)+1
        self.setup_vectors()

    def change_dt(self, dt):
        self.simParams.dt = dt
        self.simParams.num_timesteps = int(self.simParams.tstop/self.simParams.dt)+1
        self.setup_vectors()

    def setup_vectors(self):
        if self.Vm_vec[0].size() != self.simParams.num_timesteps:
            self.Vm_vec = [v.resize(self.simParams.num_timesteps) for v in self.Vm_vec]
            self.Im_vec = [i.resize(self.simParams.num_timesteps) for i in self.Im_vec]
    
    def setup_sim(self):
        h.dt = self.simParams.dt
        h.tstop = self.simParams.tstop
        h.celsius = self.simParams.celsius
        h.finitialize(self.simParams.v_init)

    def run_sim(self):
        self.setup_sim()
        
        # cache some parameters
        elec_dist_um_y = self.simParams.elec_dist * 10**3
        rho_e = self.nerve.params.rho_e
        if type(self.nerve) is MyelinatedNerve:
            #print "MyelinatedNerve!"
            elec_node_idx = self.simParams.elec_node_idx
            num_nodes = self.nerve.params.num_nodes
            if type(self.elec) is DummyBipolarElectrode:
                #print "BipolarElec!"
                # need two r_vec
                offset1 = self.simParams.node_offset
                offset2 = self.elec.dsep*1000    # [um] sep between bipolar elec
                (r_vec1, r_vec2)=self.nerve.calcDist(elec_node_idx, elec_dist_um_y, offset1, offset2)
            else:
                offset1 = self.simParams.node_offset
                (r_vec, _) = self.nerve.calcDist(elec_node_idx, elec_dist_um_y, offset1,0)
        elif type(self.nerve) is UnmyelinatedNerve:
            elecLocX = self.simParams.elecLocX
            if type(self.elec) is DummyBipolarElectrode:
                # need two r_vec
                offset1 = self.simParams.node_offset
                offset2 = self.elec.dsep*1000
                (r_vec1, r_vec2) = self.nerve.calcDist(elecLocX, elec_dist_um_y, offset1, offset2)
            else:
                offset1 = self.simParams.node_offset
                (r_vec, _) = self.nerve.calcDist(elecLocX, elec_dist_um_y, offset1, 0)

        '''
        If dummy_stim.i is in nA, then
        I/(4*pi*sigma*r) = (I*rho_e)/(4*pi*r) = ([nA][ohm-cm])/([um]) 
                    -> multiply by (1/10^6)/(1/10^4)=10^-2 to get mV
        '''
        while(h.t < h.tstop):
            # apply extracellular potential
            amp = self.elec.dummy_stim.i
            if type(self.nerve) is MyelinatedNerve:
                if type(self.elec) is DummyBipolarElectrode:
                    for i in range(num_nodes):
                        # add e-fields from both electrodes
                        e1 = (10**-2)*rho_e*amp/(4*pi*r_vec1[i])
                        e2 = (10**-2)*rho_e*(-amp)/(4*pi*r_vec2[i])
                        self.nerve.axons[i].e_extracellular = e1+e2
                else:
                    for i in range(num_nodes):
                        self.nerve.axons[i].e_extracellular = (10**-2)*rho_e*amp/(4*pi*r_vec[i])
                #if h.t-6.0 >= 0.1:
                #    import pdb; pdb.set_trace()

            elif type(self.nerve) is UnmyelinatedNerve:
                i = 0
                for seg in self.nerve.axon:
                    if type(self.elec) is DummyBipolarElectrode:
                        # add both e-fields
                        e1 = (10**-2)*rho_e*amp/(4*pi*r_vec1[i])
                        e2 = (10**-2)*rho_e*(-amp)/(4*pi*r_vec2[i])

                    else:
                        seg.e_extracellular = (10**-2)*rho_e*amp/(4*pi*r_vec[i])
                    i = i+1
            h.fadvance()

#--------------- Simulation control utilities ----------
def AP_status(Vm_vec_py):
    """
    Action potential detection:
    First 20 nodes have maximum Vm > +8mV
    Maximum at node 19 to 21 less than 13.8mV
    """
    if np.any(np.isnan(Vm_vec_py)):
        # Too high causes nan
        return 1
    try:
        AP = sum(np.max(Vm_vec_py[0:21],1)>=8)
        #print "     %d nodes higher than 8mV" % (AP)
        AP = (AP>=15)
    except:
        import pdb; pdb.set_trace()
        print "here"
    #print "max(Vm_vec_py[14:21])=", np.max(Vm_vec_py[14:21])
    tooHigh = (np.abs(np.max(Vm_vec_py[14:21]) - 13) > 1.0)

    if not AP:
        #print "     No AP, too low!"
        return -1
    elif AP and not tooHigh:
        #import pdb; pdb.set_trace()
        return 0
    elif AP and tooHigh:
        #print "     Got AP, too high"
        return 1

def AP_msgFn(mySim, myamp, niter):
    fiber_diam = mySim.nerve.params.fiber_diam
    elec_dist = mySim.simParams.elec_dist
    node_offset = mySim.simParams.node_offset/mySim.nerve.params.internodal_len
    if type(mySim.elec) is DummyBipolarElectrode:
        elec_sep = mySim.elec.dsep
    else:
        elec_sep = 0
    dur = mySim.elec.dummy_stim.dur
    if ~np.isnan(myamp):
        print "D=%.2fum, r=%.2fmm, dur=%.2fms, offset=%.2f, spacing=%.2fmm: %.2fnA or %.2fmA see action potential in %d iterations" % \
        (fiber_diam, elec_dist, dur, node_offset, elec_sep, myamp, myamp/10**6, niter)
    else:
        print "D=%.2fum, r=%.2fmm, dur=%.2fms, offset=%.2f, spacing=%.2fmm: failed to get action potential" % (fiber_diam, elec_dist, dur, node_offset, elec_sep)

def runSimAgain(mySim, statusFn):
    mySim.run_sim()
    Vm_vec_py = np.array(mySim.Vm_vec)
    status = statusFn(Vm_vec_py)
    return (Vm_vec_py, status)

def strength_duration(mySim, durations, lo_lim, hi_lim, evalFn, msgFn):
    """
    Vary strength and duration, for fixed diameter
    """
    threshCur = []
    for dur in durations:
        lo_amp = lo_lim
        hi_amp = hi_lim

        if len(threshCur) > 0:
            last_not_nan = np.where(~np.isnan(threshCur))[0]
            if len(last_not_nan)>0:
                hi_amp = np.min([threshCur[last_not_nan[-1]]*1.05, -0.5e-3/1e-9])
                print "hi_amp from cached: %.2fnA" % (hi_amp)

        mySim.elec.setDur(dur)
        mySim.change_tstop(np.min([10, 2+dur]))
        myamp = np.mean((lo_amp,hi_amp))
        mySim.elec.setAmp(myamp)
        Vm_vec_py, status = runSimAgain(mySim, evalFn)
       
        niter = 0
        while status != 0:
            if status == 1: # current too high
                hi_amp = myamp
                #if niter>30:
                #    print "Early break, inexact threshold; myamp=%2.fnA" % (myamp)
                #    break
            elif status == -1: # current too low
                lo_amp = myamp
                #if niter>30:
                #    print "Early break, can't reach threshold: loamp=%2.fnA" % (myamp)
                #    myamp = np.nan
                #    break
            if np.abs(hi_amp) <= np.abs(lo_amp):
               print "hi_amp (%2.fnA) smaller than lo_amp (%.2fnA) now!" % (hi_amp, lo_amp)
               myamp = np.nan
               break

            myamp = np.mean((lo_amp, hi_amp))
            mySim.elec.setAmp(myamp)
            Vm_vec_py, status = runSimAgain(mySim, evalFn)
            niter = niter+1
            if niter > 10:
                print "Fail to converge after %d iterations. loamp=%2.fnA, hiamp=%.2fnA" % (niter, lo_amp, hi_amp)
                myamp = np.nan
                break
        
        msgFn(mySim, myamp, niter)
        if np.isnan(myamp) and len(np.where(~np.isnan(threshCur))[0])>0:
            print "Using the previous duration's threshCur!"
            threshCur.append(threshCur[-1])
        else:
            threshCur.append(myamp)
    return threshCur

def strength_diam(mySim, diameters, lo_lim, hi_lim, evalFn, msgFn):
    """
    Vary strength and diameter, for fixed duration
    lo_lim and hi_lim are current limit, in [nA]
    """
    threshCur = []
    cur_range = np.abs(hi_lim-lo_lim)
    tol_range = 0.001*cur_range
    for diam in diameters:
        print " diameter=%.2fum" % (diam)
        mySim.nerve.set_fiber_diam(diam)
        
        lo_amp = lo_lim
        hi_amp = hi_lim
        if len(threshCur) > 0 and not np.isnan(threshCur[-1]):
            hi_amp = threshCur[-1]
        
        myamp = np.mean((lo_amp, hi_amp))
        print " setting amp=%.2f nA" % (myamp)
        mySim.elec.setAmp(myamp)
        
        Vm_vec_py, status = runSimAgain(mySim, evalFn)

        while status != 0:
            if status == 1: # current too high
                hi_amp = myamp
            elif status == -1: # current too low
                lo_amp = myamp
            #elif np.isnan(status):
            #    print "Settings cannot converge!"
            #    break

            myamp = np.mean((lo_amp, hi_amp))
            if (np.abs(hi_lim-myamp)<tol_range) or \
               (np.abs(myamp-lo_lim)<tol_range):
                print " Current limit reached, not converging!"
                break
            elif (np.abs(hi_amp-lo_amp) < tol_range) and status == 1:
                print " Stopping optimization, setting status=0"
                status = 0
            else:
                #print " setting amp=%.2f nA" % (myamp)
                mySim.elec.setAmp(myamp)
                Vm_vec_py, status = runSimAgain(mySim, evalFn)

        if status == 0:
            msgFn(mySim)
            threshCur.append(myamp)
        else:
            threshCur.append(np.nan)

    return threshCur

