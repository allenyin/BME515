from PNS import *

# Initialize sweeney nerve
initBunch = Bunch()
initBunch.fiber_diam = 12.0 # [um]
initBunch.scale_diam = 0.7
initBunch.num_nodes = 51
initBunch.num_seg = 1
initBunch.rho_a = 100.0      # [ohm-cm]
initBunch.rho_e = 500.0      # [ohm-cm]
initBunch.node_cm = 1.0     # [uF/cm^2]
initBunch.node_len = 1.5    # [um]
initBunch.scale_len = 100.0
nerve = MyelinatedNerve(initBunch, sweeneyCh)  

# Initialize passive nerve with passive characteristics
initBunch2 = Bunch()
initBunch2.fiber_diam = 12.0 # [um]
initBunch2.scale_diam = 0.7
initBunch2.num_nodes = 51
initBunch2.num_seg = 1
initBunch2.rho_a = 100.0      # [ohm-cm]
initBunch2.rho_e = 500       # [ohm-cm]
initBunch2.node_cm = 1.0     # [uF/cm^2]
initBunch2.node_len = 1.5    # [um]
initBunch2.scale_len = 100.0
nerve2 = MyelinatedNerve(initBunch2, passiveSweeneyCh)  

# Initialize stimulating electrode
myamp = -0.115e-3/1e-9
mydur = 0.2 # ms
mydel = 5
elec = DummyElectrode(mydel, myamp, mydur)

# Initialize simulation
simBunch = Bunch()
simBunch.dt = 0.01     # [ms], used to be 0.005
simBunch.tstop = 10.0   # [ms]
simBunch.v_init = -80   # [mV]
simBunch.celsius = 37   # [deg C]
simBunch.elec_node_idx = initBunch.num_nodes/2
simBunch.elec_dist = 1  # [mm]
sim1 = StimSim(nerve, elec, simBunch)   # active Sweeney
sim2 = StimSim(nerve2, elec, simBunch)  # passive Sweeney

def testPlot(sim):
    """
    Run sim and plot time course along fiber during stimulation,
    given a simulation object
    """
    sim.run_sim()
    Vm_vec_py = np.array(sim.Vm_vec)
    plt.figure()
    for i in np.arange(1000, 2000, 10):
        plt.plot(np.arange(1,52), Vm_vec_py[:,i], 'r-')
    plt.show(block=False)
    
def AP_status(Vm_vec_py):
    """
    Action potential detection:
    50% of the nodes' maximum Vm > +8mV
    The global maximum needs to be less than 13.8mV
    """
    try:
        AP = ((sum(np.max(Vm_vec_py, 1) >= 8)/(Vm_vec_py.shape[0]*1.0)) >= 0.5)
    except:
        import pdb; pdb.set_trace()
        print "here"
    print "max(Vm_vec_py) is ", np.max(Vm_vec_py)
    tooHigh = (np.max(Vm_vec_py) - 13 > 0.8)

    if not AP:
        return -1
    elif AP and not tooHigh:
        return 0
    else:
        return 1

def AP_msgFn(mySim):
    fiber_diam = mySim.nerve.params.fiber_diam
    elec_dist = mySim.simParams.elec_dist
    dur = mySim.elec.dummy_stim.dur
    myamp = mySim.elec.dummy_stim.amp
    print "D=%.2fum, r=%.2fmm, dur=%.2fms: %.2fnA or %.2fmA see action potential" % \
          (fiber_diam, elec_dist, dur, myamp, myamp/10**6)

def passive_status(Vm_vec_py):
    peak = np.max(Vm_vec_py)
    if (np.abs(peak - (sim2.simParams.v_init+20))) <= 0.1:
        return 0
    elif peak > (sim2.simParams.v_init+20):
        return 1
    else:
        return -1

def passive_msgFn(mySim):
    fiber_diam = mySim.nerve.params.fiber_diam
    elec_dist = mySim.simParams.elec_dist
    dur = mySim.elec.dummy_stim.dur
    myamp = mySim.elec.dummy_stim.amp
    dVm = np.max(np.array(mySim.Vm_vec)) - sim2.simParams.v_init
    print "D=%.2fum, r=%.2fmm, dur=%.2fms: %.2fnA or %.2fmA see dVm=%.2fmV" % \
          (fiber_diam, elec_dist, dur, myamp, myamp/10**6, dVm)

def runSimAgain(mySim, statusFn):
    mySim.run_sim()
    Vm_vec_py = np.array(mySim.Vm_vec)
    status = statusFn(Vm_vec_py)
    return (Vm_vec_py, status)

def strength_duration(mySim, durations, lo_lim, hi_lim, evalFn, msgFn):
    threshCur = []
    for dur in durations:
        lo_amp = lo_lim
        hi_amp = hi_lim

        if len(threshCur) > 0:
            hi_amp = threshCur[-1]
        mySim.elec.setDur(dur)
        myamp = np.mean((lo_amp,hi_amp))
        mySim.elec.setAmp(myamp)
        Vm_vec_py, status = runSimAgain(mySim, evalFn)
       
        while status != 0:
            if status == 1: # current too high
                hi_amp = myamp
            elif status == -1: # current too low
                lo_amp = myamp
            myamp = np.mean((lo_amp, hi_amp))
            mySim.elec.setAmp(myamp)
            Vm_vec_py, status = runSimAgain(mySim, evalFn)
        msgFn(mySim)
        threshCur.append(myamp)
    return threshCur

#------------------ Problem solving starts -------------------
durations = np.array([10e-6, 20e-6, 50e-6, \
                      100e-6, 200e-6, 500e-6, \
                      1e-3, 2e-3, 5e-3, 10e-3])     # in seconds
durations = durations/1e-3  # convert to ms
def p1_1():
    """
    Strength-duration curve for Active Sweeney nerve with cathodic stim
    """
    lo_lim = -0.05e-3/1e-9      # cathodic current low limit [nA]
    hi_lim = -3e-3/1e-9         # cathodic current high limit [nA]

    print "\nStrength-duration curve for Active Sweeney nerve with cathodic stim"
    threshCur = strength_duration(sim1, durations, lo_lim, hi_lim, AP_status, AP_msgFn)
    plt.figure()
    plt.semilogx(durations, np.abs(threshCur)/10**6, '*-')
    plt.xlabel('Duration (ms)')
    plt.ylabel('Threshold current magnitude (mA)')
    plt.title("Active Sweeney Strength-Duration, Cathodic Stim")
    plt.show(block=False)
    plt.savefig('sweeneyActive_strengthDuration_cathodic.png', bbox_inches='tight')
    return threshCur

def p1_2():
    """
    Strength-duration curve for Passive Sweeney nerve with cathodic stim
    """
    lo_lim = -0.01e-3/1e-9      # cathodic current low limit [nA]
    hi_lim = -3e-3/1e-9         # cathodic current high limit [nA]

    print "\nStrength-duration curve for Passive Sweeney nerve with cathodic stim"
    threshCur = strength_duration(sim2, durations, lo_lim, hi_lim, passive_status, passive_msgFn)
    plt.figure()
    plt.semilogx(durations, np.abs(threshCur)/10**6, '*-')
    plt.xlabel('Duration (ms)')
    plt.ylabel('Threshold current magnitude (mA)')
    plt.title("Passive Sweeney Strength-Duration, Cathodic Stim")
    plt.show(block=False)
    plt.savefig('sweeneyPassive_strengthDuration_cathodic.png', bbox_inches='tight')
    return threshCur

def p1_3_1():
    """
    Strength-duration curve for Active Sweeney nerve with anodic stim
    """
    lo_lim = 0.05e-3/1e-9      # anodic current low limit [nA]
    hi_lim = 3e-3/1e-9         # anodic current high limit [nA]

    print "\nStrength-duration curve for Active Sweeney nerve with anodic stim"
    threshCur = strength_duration(sim1, durations, lo_lim, hi_lim, AP_status, AP_msgFn)
    plt.figure()
    plt.semilogx(durations, np.abs(threshCur)/10**6, '*-')
    plt.xlabel('Duration (ms)')
    plt.ylabel('Threshold current magnitude (mA)')
    plt.title("Active Sweeney Strength-Duration, Anodic Stim")
    plt.show(block=False)
    plt.savefig('sweeneyActive_strengthDuration_anodic.png', bbox_inches='tight')
    return threshCur

def p1_3_2(threshCurCat, threshCurAn):
    """
    Compare the action potential invoked with 200us stimulation in
    cathodic vs. anodic
    """
    amp_cat = threshCurCat[4]
    amp_an = threshCurAn[4]
    #sim1.elec.setDur(200e-6/1e-3)   # 200us stimulation
    sim1.elec.setDur(1)
    
    # cathodic stimulation
    sim1.elec.setAmp(amp_cat)
    sim1.run_sim()
    Vm_cat = np.array(sim1.Vm_vec)
    
    # anodic stimulation
    sim1.elec.setAmp(amp_an)
    sim1.run_sim()
    Vm_an = np.array(sim1.Vm_vec)

    f, axarr = plt.subplots(2, sharex=True)
    """
    stimulation starts at 5ms and ends at 10ms. This is sample 1000 to 2000
    The stimulation pulse is 200us=0.2ms, this takes 0.2/0.005=40 samples.
    Let's plot from sample 1000 to 1000+2*40=1080.
    This corresponds to time from 5ms to 5.4ms
    """
    for i in np.arange(1000, 1081, 4):
        alpha = 1-(i-1000)/81.0
        axarr[0].plot(np.arange(1,52), Vm_cat[:,i], 'b-', alpha=alpha)
        axarr[1].plot(np.arange(1,52), Vm_an[:,i], 'r-', alpha=alpha)
    axarr[0].set_title('Cathodic stimulation, 5ms to 5.4ms')
    axarr[0].set_xlim([0, 51])
    axarr[0].set_ylabel('Vm (mV)')

    axarr[1].set_title('Anodic stimulation, 5ms to 5.4ms')
    axarr[1].set_xlim([0, 51])
    axarr[1].set_xlabel('Node number')
    axarr[1].set_ylabel('Vm (mV)')
    plt.show(block=False)
    plt.savefig('sweeneyActive_strengthDuration_catVsAn.png', bbox_inches='tight')
    return (Vm_cat, Vm_an)

def p1_4(Vm_cat, Vm_an, mySim):
    """
    Calculate conduction velocity from Vm of the simulation
    """
    f, axarr = plt.subplots(2, sharex=True)
    startNode = 10
    endNode = 50    # not 51 to avoid edge effects 
    catMaxIdx = []
    anMaxIdx = []
    for i in np.arange(startNode, endNode, 2):
        alpha = 1-(i-startNode)/(endNode-startNode*1.0)
        axarr[0].plot(mySim.t_vec, Vm_cat[i,:], 'b-', alpha=alpha)
        catMaxIdx.append(np.argmax(Vm_cat[i,:]))
        
        axarr[1].plot(mySim.t_vec, Vm_an[i,:], 'r-', alpha=alpha)
        anMaxIdx.append(np.argmax(Vm_an[i,:]))
    
    axarr[0].set_title('Cathodic stimulation AP, node %d to %d' % (startNode, endNode))
    axarr[0].set_ylabel('Vm (mV)')
    axarr[1].set_title('Anodic stimulation AP, node %d to %d' % (startNode, endNode))
    axarr[1].set_xlabel('Time (mS)')
    axarr[1].set_ylabel('Vm (mV)')
    plt.show(block=False)
    #plt.savefig('sweeneyActive_strengthDuration_catVsAn_time.png', bbox_inches='tight')

    # calculate conduction velocity
    dist = mySim.nerve.params.internodal_len * 2    # [um]
    avgCatTime = np.mean(np.diff(catMaxIdx)) * mySim.simParams.dt   # [mS]
    avgAnTime = np.mean(np.diff(anMaxIdx)) * mySim.simParams.dt     # [mS]

    print "Cathodic conduction velocity is %.2fm/s" % ((dist/1e6)/(avgCatTime/1e3))
    print "Anodic conduction velocity is %.2fm/s" % ((dist/1e6)/(avgAnTime/1e3))

def p1_5(Ith, durations):
    """
    Fit the Weiss form of strength-duration curve:
    Ith(PW) = Irh(1+Tch/PW) with linear regression

    Solving by doing: Ith = Irh + (Irh*Tch)/PW
                          = Irh + (Irh*Tch)*x,  x=1/PW

    Then linear regression gives Irh and Irh*Tch
    """
    n = len(durations)
    Ith = np.array(Ith)
    Ith = Ith.reshape((n,1))
    A = np.ones((n,2))
    A[:,1] = 1/durations
    coeffs,resid,rank,s = np.linalg.lstsq(A, Ith)
    Irh = coeffs[0]     # in nA
    Tch = coeffs[1]/Irh # in ms
    print "\nFit strength-duration curve for Active Sweeney nerve with Cathodic stim"
    print "Irh=%.2fmA, Tch=%.2fmS" % (Irh/10**6, Tch)

    weissFn = lambda x: Irh*(1+Tch/x)       # gives value in nA
    # calculate percent residuals
    Ith_pred = [weissFn(i) for i in durations]
    perc_residuals = [np.abs((Ith_pred[i]-Ith[i])/Ith[i]) for i in range(len(Ith))]
    print "Average residual of the fit is %.2f%%" % (np.mean(perc_residuals)*100)

    # make the fit for plotting
    fitDur = np.logspace(-2, 1, num=100)   # logspace from 0.01ms to 10ms
    fit = [weissFn(i) for i in fitDur]
    fit = np.concatenate(fit, axis=1)
    Ich_ab = np.abs(weissFn(Tch))/10**6
    Irh_ab = np.abs(Irh)/10**6

    plt.figure()
    l1, = plt.semilogx(durations, np.abs(Ith)/10**6, 'b*-')
    l2, = plt.semilogx(fitDur, np.abs(fit)/10**6, 'r-')

    # drawing lines go from (x1,x2) to (y1,y2)..
    l3, = plt.semilogx((fitDur[0], fitDur[-1]), (Irh_ab, Irh_ab), 'r--')    # rheobase
    l4, = plt.semilogx((Tch, Tch), (0, Ich_ab), 'b--')            # chronaxie vertical
    l5, = plt.semilogx((fitDur[0], Tch), (Ich_ab, Ich_ab), 'b--') # chronaxie horizontal

    plt.xlabel('Duration (ms)')
    plt.ylabel('Threshold current magnitude (mA)')
    plt.title("Active Sweeney Strength-Duration, Cathodic Stim")
    plt.legend([l1,l2,l3,l4], ['Simulation', 'Weiss Fit', 'Rheobase', 'Chronaxie'])
    plt.show(block=False)
    plt.savefig('sweeneyActive_cathodic_simVsFit.png', bbox_inches='tight')
    return (Irh, Tch)

    

threshCur1 = p1_1()
threshCur2 = p1_2()
threshCur3 = p1_3_1()
#(Vm_cat, Vm_an) = p1_3_2(threshCur1, threshCur3)
#p1_4(Vm_cat, Vm_an, sim1)
#Irh, Tch = p1_5(threshCur1, durations)
