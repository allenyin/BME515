from PNS import *
import numpy as np
import pdb

# Run simulation one fiber at a time

rate = 'baseRate'   # or 'activeRate'
#rate = 'activeRate'
saveDone = True 
calcSNR = True

#--------------------------
pop1 = [15, 2.5]    # alpha, beta motor
pop2 = [4, 1.5]     # pain fibers

np.random.seed(13)  # set seed so we get fixed results
bundle_diam = np.concatenate([np.random.normal(15, 2.5, 50), np.random.normal(4, 1.5, 50)])
base_freq = np.random.randint(10, 20, 100)
active_freq = np.random.randint(30, 60, 100)

# Initialize sweeney nerve parameters
initBunch = Bunch()
initBunch.scale_diam = 0.6
initBunch.num_nodes = 51    # number of nodes the longest fiber will have
initBunch.fiber_diam = 12   # default...will change
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

We want at 51 nodes for the longest nerve
'''
max_internodal_len = initBunch.scale_len * bundle_diam[np.argmax(bundle_diam)]+0.5*initBunch.node_len
max_nerve_len = max_internodal_len * initBunch.num_nodes

# Get the number of nodes needed for each fiber
list_num_nodes = []
for i in np.arange(len(bundle_diam)):
    cur_internodal_len = initBunch.scale_len*bundle_diam[i] + 0.5*initBunch.node_len
    list_num_nodes.append(int(max_nerve_len/cur_internodal_len)//2*2+1)

'''
Simulation setup

Lowest frequency is 10Hz, that's 100ms/pulse. 
Want at least 3 pulses, So simulation will take at least 300ms

Let dt=0.015ms, so the integration doesn't fuck up 
'''
tstop = 250 # [ms]
dt = 0.015    # [ms]
v_init = -80.0 # [mV]
celsius = 37.0 # [deg C]
num_timesteps = int(tstop/dt)+1

#----------------------------------------------------------------------------------
# Do simulation if we haven't collected data yet
#
# Run simulation one fiber at a time, then save data
#----------------------------------------------------------------------------------
if not saveDone:
    print 'Collecting data for %s' % (rate)
    
    for i in range(len(bundle_diam)):
        initBunch.fiber_diam = bundle_diam[i]
        initBunch.num_nodes = list_num_nodes[i]
        nerve = MyelinatedNerve(initBunch, sweeneyCh)
        stim = h.Ipulse2(nerve.axons[0](0.5))

        stim.delay = 2    # [ms]
        stim.dur = 0.1    # [ms]
        stim.amp = 3.5    # [nA], enough for all to fire
        stim.num = 20     # In 300ms, w/ f_max=60, we need max 18 pulses                

        if rate is 'baseRate':
            stim.per = 1000.0/base_freq[i]
        elif rate is 'activeRate':
            stim.per = 1000.0/active_freq[i]
        
        Vm_vec = []
        Im_vec = []
        for j in range(nerve.params.num_nodes):
            Vm_vec.append(h.Vector(num_timesteps,0))
            Vm_vec[j].record(nerve.axons[j](0.5)._ref_v)

            Im_vec.append(h.Vector(num_timesteps,0))
            Im_vec[j].record(nerve.axons[j](0.5)._ref_i_membrane)
        t_vec = h.Vector()
        t_vec.record(h._ref_t)

        # start sim
        h.dt = dt
        h.tstop = tstop
        h.celsius = celsius
        h.finitialize(v_init)
        while (h.t < h.tstop):
            h.fadvance()

        # Save the data
        Vm_vec_py = np.array(Vm_vec)
        Im_vec_py = np.array(Im_vec)
        np.save('./sim_data/Vm_vec_py_fiber%d_%s.npy' % (i, rate), Vm_vec_py)
        np.save('./sim_data/Im_vec_py_fiber%d_%s.npy' % (i, rate), Im_vec_py)
    
    # Save time just once
    t_vec_py = np.array(t_vec)
    np.save('./sim_data/t_vec_py_%s.npy' % (rate), t_vec_py)

elif saveDone and calcSNR:
    print 'Calculate SNR'

    K = [5.0, 10.0, 15.0, 20.0]
    V1_big_base = []
    V2_big_base = []
    V3_big_base = []
    V1_small_base = []
    V2_small_base = []
    V3_small_base = []

    V1_big_active = []
    V2_big_active = []
    V3_big_active = []
    V1_small_active = []
    V2_small_active = []
    V3_small_active = []

    SNR_small = []
    SNR_big = []
    for kIdx in range(len(K)):
        npzfile = np.load('Vext_k=%d_baseRate.npz' % (int(K[kIdx])))
        V1_big_base.append(npzfile['arr_0'])
        V2_big_base.append(npzfile['arr_1'])
        V3_big_base.append(npzfile['arr_2'])
        V1_small_base.append(npzfile['arr_0'])
        V2_small_base.append(npzfile['arr_1'])
        V3_small_base.append(npzfile['arr_2'])
        
        npzfile = np.load('Vext_k=%d_activeRate.npz' % (int(K[kIdx])))
        V1_big_active.append(npzfile['arr_0'])
        V2_big_active.append(npzfile['arr_1'])
        V3_big_active.append(npzfile['arr_2'])
        V1_small_active.append(npzfile['arr_0'])
        V2_small_active.append(npzfile['arr_1'])
        V3_small_active.append(npzfile['arr_2'])

        V_big_base_ref = (V1_big_base[kIdx] + V3_big_base[kIdx])/2
        V_big_base_signal = V2_big_base[kIdx] - V_big_base_ref

        V_big_active_ref = (V1_big_active[kIdx] + V3_big_active[kIdx])/2
        V_big_active_signal = V2_big_active[kIdx] - V_big_active_ref

        V_small_base_ref = (V1_small_base[kIdx] + V3_small_base[kIdx])/2
        V_small_base_signal = V2_small_base[kIdx] - V_small_base_ref

        V_small_active_ref = (V1_small_active[kIdx] + V3_small_active[kIdx])/2
        V_small_active_signal = V2_small_active[kIdx] - V_small_active_ref

        # SNR = 20*log(variance(active_signal)/(mse(active_signal-base_signal)))
        mse_big = np.mean((V_big_active_signal - V_big_base_signal)**2)
        var_big = np.var(V_big_active_signal)
        SNR_big.append(20*np.log10(var_big/mse_big))

        mse_small = np.mean((V_small_active_signal - V_small_base_signal)**2)
        var_small = np.var(V_small_active_signal)
        SNR_small.append(20*np.log10(var_small/mse_small))

    plt.figure()
    plt.plot(K, SNR_big, 'b*-')
    plt.plot(K, SNR_small, 'r*-')
    plt.legend(['big', 'small'])
    plt.xlabel('Electrode spacing (mm)')
    plt.ylabel('SNR (dB)')
    plt.show(block=False)
    plt.savefig('p3_SNR.png', bbox_inches='tight')

else:
    print 'Data already collected, calculating extracellular potential..'
    elec_dist = 1 # electrode 1mm away from fibers
    nNerves = len(bundle_diam)
    t_vec_py = np.load('./sim_data/t_vec_py_%s.npy' % (rate))

    K = [5.0, 10.0, 15.0, 20.0]
    V1_big = []
    V2_big = []
    V3_big = []
    V1_small = []
    V2_small = []
    V3_small = []
    
    # Initialize electrode measurment
    for kk in K:
        V1_big.append(np.zeros(len(t_vec_py)))
        V2_big.append(np.zeros(len(t_vec_py)))
        V3_big.append(np.zeros(len(t_vec_py)))
        V1_small.append(np.zeros(len(t_vec_py)))
        V2_small.append(np.zeros(len(t_vec_py)))
        V3_small.append(np.zeros(len(t_vec_py)))

    # Load fiber one-by-one and calculate contribution to electrode
    for j in range(nNerves):
        print 'Nerve %d...' % (j)
        internodal_len = initBunch.scale_len * bundle_diam[j] + 0.5*initBunch.node_len
        num_nodes = list_num_nodes[j]
        idx_center_node = (num_nodes-1)/2   # also half-len

        cur_Im_vec = np.load('./sim_data/Im_vec_py_fiber%d_%s.npy' % (j, rate))
        node_diam = initBunch.scale_diam * bundle_diam[j]
        area = pi * node_diam * initBunch.node_len
        cur_Im_vec *= area
        
        for kIdx in range(len(K)):
            kval = K[kIdx]
            e2_offset = 0              # electrode 2 aligned with center electrode
            e1_offset = -kval*(10**3)  # electrode 1 offset from center [um]
            e3_offset = kval*(10**3)   # electrode 3 offset from center [um]

            for i in range(num_nodes):
                node_offset = (i-idx_center_node)*internodal_len
                dist1 = np.sqrt((node_offset - e1_offset)**2 + elec_dist**2)
                dist2 = np.sqrt((node_offset - e2_offset)**2 + elec_dist**2)
                dist3 = np.sqrt((node_offset - e3_offset)**2 + elec_dist**2)

                for t in range(len(t_vec_py)):
                    if j >= 50:
                        V1_small[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist1)
                        V2_small[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist2)
                        V3_small[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist3)
                    else:
                        V1_big[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist1)
                        V2_big[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist2)
                        V3_big[kIdx][t] += cur_Im_vec[i, t]/(4*pi*dist3)

    # Plot the ENG
    for kIdx in range(len(K)):
        V_ref_big = (V1_big[kIdx] + V3_big[kIdx])/2
        V_ref_small = (V1_small[kIdx] + V3_small[kIdx])/2
        V_ref_all = (V1_big[kIdx] + V3_big[kIdx] + V1_small[kIdx] + V3_small[kIdx])/2

        f, (ax1,ax2,ax3) = plt.subplots(3, sharex = True)
        ax1.plot(t_vec_py, V2_big[kIdx] - V_ref_big, 'b-')
        ax2.plot(t_vec_py, V2_small[kIdx] - V_ref_small, 'r-')
        ax3.plot(t_vec_py, (V2_big[kIdx]+V2_small[kIdx])-V_ref_all, 'k')

        ax1.set_title('ENG, big fibers, K=%.2f mm, %s' % (K[kIdx], rate))
        ax2.set_title('ENG, small fibers')
        ax3.set_title('ENG, all fibers')
        ax1.set_ylabel('V (mV)')
        ax2.set_ylabel('V (mV)')
        ax3.set_ylabel('V (mV)')
        ax3.set_xlabel('Time (mS)')
        plt.show(block=False)

        plt.savefig('p3_ENG_%s_k=%d.png' % (rate, K[kIdx]), bbox_inches='tight')

        np.savez('Vext_k=%d_%s.npz' % (int(K[kIdx]), rate), \
                V1_big[kIdx], V2_big[kIdx], V3_big[kIdx], \
                V1_small[kIdx], V2_small[kIdx], V3_small[kIdx])

 
