from PNS import *
import numpy as np
import pdb

do_save = 'V'
save_Done = True

'''
2 populations of fibers:
    myelinated motor (alpha, beta): mean=15um, sd=2.5um. So 2 s.d. = 5um
    pain fibers (type-III): mean=4um, sd=1.5um. So 2 s.d. = 3um.

    Randomly select which population a fiber belongs to, then draw from
    the corresponding population
'''
pop1 = [15, 2.5]    # alpha, beta motor
pop2 = [4, 1.5]     # pain fibers

np.random.seed(13)  # set seed so we get fixed results
bundle_diam = np.concatenate([np.random.normal(15, 2.5, 50), np.random.normal(4, 1.5, 50)])

def plot_fiber_diam():
    plt.figure()
    plt.hist(bundle_diam, range=(0,25))
    plt.xlabel('Fiber diameter (um)')
    plt.ylabel('Number of fibers')
    plt.title('Distribution of fiber diameters')
    plt.show(block=False)
    plt.savefig('fiber_diam.png', bbox_inches='tight')

# Initialize sweeney nerve bundles
initBunch = Bunch()
initBunch.scale_diam = 0.6
initBunch.num_nodes = 51
initBunch.num_seg = 1
initBunch.rho_a = 54.7      # [ohm-cm]
initBunch.rho_e = 500       # [ohm-cm]
initBunch.node_cm = 2.5     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0

'''
We want the nerves in the bundle to be the same length,
due to different fiber-diam and thus internodal length,
different nerves will will have different number of nodes.

We want at least 51 nodes for each nerve
'''
max_internodal_len = initBunch.scale_len * bundle_diam[np.argmax(bundle_diam)]+0.5*initBunch.node_len
max_nerve_len = max_internodal_len * initBunch.num_nodes

nerveBundle = []
for i in np.arange(100):
    initBunch.fiber_diam = bundle_diam[i]   # [um]
    cur_internodal_len = initBunch.scale_len * bundle_diam[i]+0.5*initBunch.node_len
    initBunch.num_nodes = int(max_nerve_len/cur_internodal_len)//2*2+1
    nerveBundle.append(MyelinatedNerve(initBunch, sweeneyCh))

# Initialize stimulating electrode
myamp = -0.5e-3/1e-9   # 0.2mA cathodic current for 0.4ms duration
mydur = 0.4 # ms
mydel = 2
elec = DummyElectrode(mydel, myamp, mydur)

# Initialize bundle simulation
#stimNodeScale = 20      # number of nodes to the left of center to stimulate
stimNodeScale = 22

simBunch = Bunch()
simBunch.dt = 0.005     # [ms]
simBunch.tstop = 8.0   # [ms]
simBunch.v_init = -80   # [mV]
simBunch.celsius = 37   # [deg C]
simBunch.stimDistOffset = -max_internodal_len * stimNodeScale # [um] stimulate at center-(this distance)
simBunch.elec_dist = 1  # [mm]

#---------------------------------------------------------------------------------
# Utility functions
#---------------------------------------------------------------------------------
def plot_nerve_vec(nerveBundle, nerveIdx, t_vec_py, measurement, m_type, tIdx):
    '''
    Stimulation starts @ 2ms, total=8ms, lasts 0.4ms, dt=0.005ms
    So plot from 1.5ms to 6ms --> samples 300 to 1200.

    Plot Vm, Im at all nodes, for time sample number [400, 420, 440, 460, 480, 500, 520] 
    x-axis should be in distance
    '''
    internodal_len = nerveBundle[nerveIdx].params.internodal_len
    num_nodes = nerveBundle[nerveIdx].params.num_nodes
    x_axis = np.arange(num_nodes) * internodal_len # [um]
   
    center_row = (measurement.shape[0]-1)/2
    half_len = (num_nodes-1)/2
    row_idx = np.arange(center_row-half_len, center_row+half_len+1)
    
    plt.figure()
    for t in tIdx:
        plt.plot(x_axis, measurement[row_idx, t, nerveIdx])

    plt.title('Nerve %d, diam=%.2fum' % (nerveIdx, nerveBundle[nerveIdx].params.fiber_diam))
    plt.ylabel(m_type)
    plt.xlabel('Dist (um)')
    plt.show(block=False)

def getExtV(e1_offset, e2_offset, e3_offset, nerveBundle, elec_dist, Im_vec_py, t_vec_py):
    '''
    Get the time course measured at extracellular electrodes by superposition
    e1_offset to e3_offset are in [um]
    elec_dist is vertical distance from electrode to nerve in [um]

    The first 50 nerves are assumed to be from different distribution than
    the second 50 nerves
    '''
    V1_big = np.zeros(len(t_vec_py))
    V2_big = np.zeros(len(t_vec_py))
    V3_big = np.zeros(len(t_vec_py))
    V1_small = np.zeros(len(t_vec_py))
    V2_small = np.zeros(len(t_vec_py))
    V3_small = np.zeros(len(t_vec_py))

    center_row = (Im_vec_py.shape[0]-1 )/2
    for j in range(len(nerveBundle)):
        internodal_len = nerveBundle[j].params.internodal_len
        num_nodes = nerveBundle[j].params.num_nodes
        idx_center_node = (num_nodes-1)/2   # also half-len
        for i in range(num_nodes):
            node_offset = (i-idx_center_node)*internodal_len
            dist1 = np.sqrt((node_offset - e1_offset)**2 + elec_dist**2)
            dist2 = np.sqrt((node_offset - e2_offset)**2 + elec_dist**2)
            dist3 = np.sqrt((node_offset - e3_offset)**2 + elec_dist**2)
            node_row = center_row-idx_center_node+i # find correct entry in array
            for t in range(len(t_vec_py)):
                if Im_vec_py[node_row,t,j] is not np.nan:
                    if j >= 50:
                        V1_small[t] = V1_small[t] + Im_vec_py[node_row, t, j]/(4*pi*dist1)
                        V2_small[t] = V2_small[t] + Im_vec_py[node_row, t, j]/(4*pi*dist2)
                        V3_small[t] = V3_small[t] + Im_vec_py[node_row, t, j]/(4*pi*dist3)
                    else:
                        V1_big[t] = V1_big[t] + Im_vec_py[node_row, t, j]/(4*pi*dist1)
                        V2_big[t] = V2_big[t] + Im_vec_py[node_row, t, j]/(4*pi*dist2)
                        V3_big[t] = V3_big[t] + Im_vec_py[node_row, t, j]/(4*pi*dist3)
    return V1_big, V2_big, V3_big, V1_small, V2_small, V3_small

def plotCNAP(V_ref, V_active, t_vec_py, K):
    plt.figure()
    plt.plot(t_vec_py, V_active-V_ref, 'b-')
    plt.title('CNAP spacing=%.2f mm' % (K))
    plt.xlabel('Time (mS)')
    plt.ylabel('Extracellular voltage (mV)')
    plt.show(block=False)

def plotCNAP_pop(V1_big, V2_big, V3_big, V1_small, V2_small, V3_small, t_vec_py, K):
    V_ref_big = (V1_big+V3_big)/2
    V_ref_small = (V1_small+V3_small)/2
    V_ref_all = (V1_big+V3_big+V1_small+V3_small)/2

    f, (ax1,ax2,ax3) = plt.subplots(3, sharex = True)
    ax1.plot(t_vec_py, V2_big-V_ref_big, 'b-')
    ax2.plot(t_vec_py, V2_small-V_ref_small, 'r-')
    ax3.plot(t_vec_py, (V2_big+V2_small)-V_ref_all, 'k')

    ax1.set_title('CNAP, big fibers, K=%.2f mm' % (K))
    ax2.set_title('CNAP, small fibers')
    ax3.set_title('CNAP, all fibers')
    ax1.set_ylabel('V (mV)')
    ax2.set_ylabel('V (mV)')
    ax3.set_ylabel('V (mV)')
    ax3.set_xlabel('Time (mS)')
    plt.show(block=False)

def plotAtElec(V1_big, V2_big, V3_big, V1_small, V2_small, V3_small, t_vec_py, K):
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(t_vec_py, V1_big, 'b-')
    ax1.plot(t_vec_py, V2_big, 'r-')
    ax1.plot(t_vec_py, V3_big, 'g-')
    ax1.set_title('Big fibers, at electrodes, K=%.2f mm' % (K))
    ax1.set_ylabel('V (mV)')

    ax2.plot(t_vec_py, V1_small, 'b-')
    ax2.plot(t_vec_py, V2_small, 'r-')
    ax2.plot(t_vec_py, V3_small, 'g-')
    ax2.set_title('Small fibers, at electrodes, K=%.2f mm' % (K))
    ax2.set_ylabel('V (mV)')
    ax2.set_xlabel('Time (mS)')
    plt.show(block=False)

#---------------------------------------------------------------------------------

#subset = [nerveBundle[0], nerveBundle[-1]]
subset = nerveBundle

if not save_Done:
    # Do simulation if no data yet
    bundleSim1 = StimSimBundle(subset, elec, simBunch, record=do_save)
    bundleSim1.run_sim()
    t_vec_py = np.array(bundleSim1.t_vec)
    Vm_vec_py,Im_vec_py = bundleSim1.getBundleVecs()

    if Vm_vec_py is not None:
        measurement = Vm_vec_py
        np.save('Vm_vec_py_cathodicStim_%d.npy' % (stimNodeScale), Vm_vec_py)
    else:
        measurement = Im_vec_py
        np.save('Im_vec_py_cathodicStim_%d.npy' % (stimNodeScale), Im_vec_py)

    np.save('t_vec_py_cathodicStim.npy', t_vec_py)

    if do_save is 'V':
        ylabel = 'Vm (mV)'
        fname = 'Vm'
    elif do_save is 'I':
        ylabel = 'Im (mV)'
        fname = 'Im'

    tIdx = [400,440,480,520,560,600,640,680]
    plot_nerve_vec(subset, 0, t_vec_py, measurement, ylabel, tIdx)
    plt.savefig('p2_%s_fiber0_%d.png' % (fname, stimNodeScale), bbox_inches='tight')
    plot_nerve_vec(subset, -1, t_vec_py, measurement, ylabel, tIdx)
    plt.savefig('p2_%s_fiber-1_%d.png' % (fname, stimNodeScale), bbox_inches='tight')

else:
    print 'Simulation finished, calculating extracellular voltage'

    # if we have already finished the simulation, calculate measured potential
    Im_vec_py = np.load('Im_vec_py_cathodicStim_%d.npy' % (stimNodeScale))
    t_vec_py = np.load('t_vec_py_cathodicStim.npy')

    '''
    Center electrode is at the center of the fibers, which are all aligned.
    Spacing K is in mm.
    '''
    K = [5.0, 10.0, 15.0, 20.0]   # recording electrode spacing in [mm]
    for kk in K:
        e2_offset = 0            # electrode 2 aligned with center electrode
        e1_offset = -kk*(10**3)  # electrode 1 offset from center [um]
        e3_offset = kk*(10**3)  # electrode 3 offset from center [um]

        V1_big,V2_big,V3_big,V1_small,V2_small,V3_small = \
            getExtV(e1_offset, e2_offset, e3_offset, subset, simBunch.elec_dist*1000, Im_vec_py, t_vec_py)
        plotCNAP_pop(V1_big,V2_big,V3_big,V1_small,V2_small,V3_small,t_vec_py,kk)
        plt.savefig('p2_CNAP_K=%d_%d.png' % (kk, stimNodeScale), bbox_inches='tight')

        plotAtElec(V1_big, V2_big, V3_big, V1_small, V2_small, V3_small, t_vec_py, kk)
        plt.savefig('p2_allElec_K=%d_%d.png' % (kk, stimNodeScale), bbox_inches='tight')




