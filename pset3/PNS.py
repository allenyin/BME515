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

# axon initialization parameter bunch for testing
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

# simulation initialization parameter bunch for testing
simBunch = Bunch()
simBunch.dt = 0.005     # [ms]
simBunch.tstop = 10.0   # [ms]
simBunch.v_init = -70.0 # [mV]
simBunch.celsius = 37.0 # [deg C]
simBunch.elec_node_idx = initBunch.num_nodes/2
simBunch.elec_dist = 1.0  # [mm]

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
        extra_axons = self.createAxons(n_new)
        extra_axons[0].connect(self.axons[-1],1,0)
        self.axons.append(extra_axons)

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

class StimSim(object):
    """
    Class to describe an extracellular stimulation simulation.
    Requires a MyelinatedNerve and a DummyElectrode object.

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

    Each instance contains:
        nerve: MyelinatedNerve object reference
        elec:  DummyElectrode object reference
        simParams: simulation parameter bunch
        Vm_vec:   Hoc vector containing Vm      (in mV)
        t_vec:    Hoc vector containing time    (in mS)
        Im_vec:   Hoc vector containing transmembrane current. (in mA/cm^2)
                  Initialized if constructor 'recordI' is true
                  
    Instead of a MyelinatedNerve, can also take a list of 
    MyelinatedNerve of different diameters. Stimulation is
    pointed at the same node number at all nerves.

    In this case, the Vm_vec and Im_vec will contain that of all MyelinatedNerves

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
        for i in range(self.nerve.params.num_nodes):
            self.Vm_vec.append(h.Vector(self.simParams.num_timesteps, 0))
            self.Vm_vec[i].record(self.nerve.axons[i](0.5)._ref_v)
            self.Im_vec.append(h.Vector(self.simParams.num_timesteps, 0))
            self.Im_vec[i].record(self.nerve.axons[i](0.5)._ref_i_membrane)
        self.t_vec = h.Vector()
        self.t_vec.record(h._ref_t)

    def getImVec(self):
        # convert the Im_vec [mA/cm^2] by multiplying area [um^2], to give nA
        # for single nerve use
        area = pi * self.nerve.params.node_diam * self.nerve.params.node_len
        print "getImVec: diam=", self.nerve.axons[0].diam
        print "getImVec: area=", area 
        Im_vec_nA = np.array(self.Im_vec)
        Im_vec_nA = Im_vec_nA*area
        return Im_vec_nA

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

    def calcDist(self, idx):
        # results in um
        elecIdx = self.simParams.elec_node_idx
        internodal_len = self.nerve.params.internodal_len
        elec_dist_um = self.simParams.elec_dist * 10**3
        return np.sqrt((np.abs(idx-elecIdx)*internodal_len)**2 + elec_dist_um**2)

    def run_sim(self):
        self.setup_sim()
        # cache some parameters
        num_nodes = self.nerve.params.num_nodes
        rho_e = self.nerve.params.rho_e

        r_vec = [self.calcDist(i) for i in range(num_nodes)] 
        '''
        If dummy_stim.i is in nA, then
        I/(4*pi*sigma*r) = (I*rho_e)/(4*pi*r) = ([nA][ohm-cm])/([um]) -> multiply by (1/10^6)/(1/10^4)=10^-2 to get mV
        '''
        while(h.t < h.tstop):
            for i in range(num_nodes):
                # apply extracellular potential
                amp = self.elec.dummy_stim.i
                self.nerve.axons[i].e_extracellular = (10**-2)*rho_e*amp/(4*pi*r_vec[i])
            h.fadvance()

class StimSimBundle(object):
    '''
    Class to describe an extracellular stimulation simulation.
    Takes in a list of MyelinatedNerve and a DummyElectrode object.

    Assume the only difference between the nerve is the diameter.
    All fibers are aligned at the center node.
    DummyElectrode is aligned to stimDistOffset to the left of the center nodes.
    Assume same distance to all nerves per pset3.

    Geometry:
                   0====0====0====0====0====0====0
                  0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0
                        ^         c 
    ^ = stimulating electrode
    c = center node
    
    simulation initialization parameter bunch example
        simBunch = Bunch()
        simBunch.dt = 0.005     # [ms]
        simBunch.tstop = 10.0   # [ms]
        simBunch.v_init = -70.0 # [mV]
        simBunch.celsius = 37.0 # [deg C]
        simBunch.stimDistOffset = 10 # [um] 
        simBunch.elec_dist = 1  # [mm]

    Each instance contains:
        bundle: List of MyelinatedNerve object reference
        elec:  DummyElectrode object reference
        simParams: simulation parameter bunch
        Vm_vec:   Hoc vector containing Vm      (in mV)
        t_vec:    Hoc vector containing time    (in mS)
        Im_vec:   Hoc vector containing transmembrane current. (in mA/cm^2)
                  Initialized if constructor 'recordI' is true
                  
    Vm_vec and Im_vec are lists (individual nerves) of lists (Vm, Im vector at each node)
    '''
    def __init__(self, bundle, elec, simBunch, record='all'):
        assert type(bundle[0]) is MyelinatedNerve
        self.bundle = bundle 
        self.elec = elec
        self.simParams = simBunch
        num_timesteps = int(simBunch.tstop/simBunch.dt)+1
        self.simParams.num_timesteps = num_timesteps

        # Record only one at a time, because memory...
        recordV = False
        recordI = False
        if record is 'all':
            recordV = True
            recordI = True
        elif record is 'V':
            recordV = True
            self.Im_vec = None
        elif record is 'I':
            recordI = True
            self.Vm_vec = None
        self.create_bundle_vectors(recordV, recordI)

    def create_bundle_vectors(self, recordV, recordI):
        if recordV:
            self.Vm_vec = []
        if recordI:
            self.Im_vec = []
        
        for j in range(len(self.bundle)):
            if recordV:
                self.Vm_vec.append([])
            if recordI:
                self.Im_vec.append([])

            for i in range(self.bundle[j].params.num_nodes):
                if recordV:
                    self.Vm_vec[j].append(h.Vector(self.simParams.num_timesteps, 0))
                    self.Vm_vec[j][i].record(self.bundle[j].axons[i](0.5)._ref_v)
                if recordI:
                    self.Im_vec[j].append(h.Vector(self.simParams.num_timesteps, 0))
                    self.Im_vec[j][i].record(self.bundle[j].axons[i](0.5)._ref_i_membrane)
        self.t_vec = h.Vector()
        self.t_vec.record(h._ref_t)

    def getBundleVecs(self):
        '''
        Convert and return Vm_vec and Im_vec as a numpy array
        Vm_vec_py and Im_vec_py will be a NxTxM array, where
            N=number of maximum nodes possible
            T=number of time points,
            M=number of nerves in the bundle

        The center nodes are aligned in the array. Since the nerves
        have different number of nodes, nan will be used as place holders
        
        Im_vec [mA/cm^2] entries will be multiplied by area [um^2], to give nA
        in result.
        '''
        num_nodes_list = [nerve.params.num_nodes for nerve in self.bundle]
        max_nodes = max(num_nodes_list)
        center_row = (max_nodes-1)/2

        if self.Vm_vec is not None:
            Vm_vec_py = np.empty((max_nodes, self.simParams.num_timesteps,len(self.bundle)))
            Vm_vec_py[:] = np.NAN
        else:
            Vm_vec_py = None

        if self.Im_vec is not None: 
            Im_vec_py = np.empty((max_nodes, self.simParams.num_timesteps,len(self.bundle)))
            Im_vec_py[:] = np.NAN
        else:
            Im_vec_py = None

        for j in range(len(self.bundle)):
            half_len = (num_nodes_list[j]-1)/2
            row_idx = np.arange((center_row-half_len), (center_row+half_len+1))

            if self.Vm_vec is not None:
                Vm_vec_py[row_idx,:,j] = np.array(self.Vm_vec[j])
            
            if self.Im_vec is not None:
                area = pi * self.bundle[j].params.node_diam * self.bundle[j].params.node_len
                Im_vec_py[row_idx,:,j] = np.array(self.Im_vec[j])
                Im_vec_py[row_idx,:,j] = Im_vec_py[row_idx,:,j] * area
        
        return (Vm_vec_py, Im_vec_py)

    def change_tstop(self, tstop):
        self.simParams.tstop = tstop
        self.simParams.num_timesteps = int(self.simParams.tstop/self.simParams.dt)+1
        self.setup_vectors()

    def change_dt(self, dt):
        self.simParams.dt = dt
        self.simParams.num_timesteps = int(self.simParams.tstop/self.simParams.dt)+1
        self.setup_vectors()

    def setup_vectors(self):
        # resize the recording vectors
        if self.Vm_vec[0][0].size() != self.simParams.num_timesteps:
            for j in range(len(self.bundle)):
                for i in range(self.bundle[j].params.num_nodes):
                    self.Vm_vec[j][i].resize(self.simParams.num_timesteps)
                    self.Im_vec[j][i].resize(self.simParams.num_timesteps)
    
    def setup_sim(self):
        h.dt = self.simParams.dt
        h.tstop = self.simParams.tstop
        h.celsius = self.simParams.celsius
        h.finitialize(self.simParams.v_init)

    def calcDist(self, nerve, idx):
        # results in um
        node_offset = (idx - (nerve.params.num_nodes-1)/2)*nerve.params.internodal_len
        dist_from_elec = np.abs(node_offset - self.simParams.stimDistOffset)    # [um]
        elec_dist_um = self.simParams.elec_dist * 10**3
        return np.sqrt(dist_from_elec**2 + elec_dist_um**2)

    def run_sim(self):
        self.setup_sim()
        
        '''
        Cache some parameters
        This simulation assumes electrode is same distance away from
        all nerves. All nerves have same number of nodes.
        '''
        num_nodes = [nerve.params.num_nodes for nerve in self.bundle]
        rho_e = [nerve.params.rho_e for nerve in self.bundle]

        '''
        Calculate r_vec -- distance from stim to each ndoes on all nerves
        r_vec is a list, each entry corresponds to a nerve.
              Each of those entries is a list of distance from node to electrode
        '''
        r_vec = []
        for j in range(len(self.bundle)):
            r_vec.append([])
            r_vec[j] = [self.calcDist(self.bundle[j], i) for i in range(num_nodes[j])]

        '''
        If dummy_stim.i is in nA, then
        I/(4*pi*sigma*r) = (I*rho_e)/(4*pi*r) = ([nA][ohm-cm])/([um]) -> multiply by (1/10^6)/(1/10^4)=10^-2 to get mV
        '''
        while(h.t < h.tstop):
            amp = self.elec.dummy_stim.i
            for j in range(len(self.bundle)):
                for i in range(num_nodes[j]):
                    self.bundle[j].axons[i].e_extracellular = \
                      (10**-2)*rho_e[j]*amp/(4*pi*r_vec[j][i])
            h.fadvance()
