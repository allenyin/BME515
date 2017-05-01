import numpy as np
from matplotlib import pyplot as plt

durations = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.2, 1.5]  # PW [ms]
elec_dist = [0.5, 0.75, 1.0, 1.25, 1.5]                     # vertical dist [mm]
node_offset = [0, 0.25, 0.5]                                # frac of internodal len
dsep = [1, 5, 10, 15, 20]                                   # lateral elec sep [mm]

#turquoise(blue), lavender(gray), green, red, mustard(yellow)
edgecolors = ['#06c2ac','#c79f3f','#06470c','#650021', '#ceb301']
facecolors = ['#75bbfd','#929591','#96f97b','#c20078', '#ffff14']

fiber_diam = [1, 3, 5,7,9]
all_data = []
for i in np.arange(len(fiber_diam)):
    all_data.append(np.load('threshCur_bipolarStim_%dum.npy' % (fiber_diam[i])))
all_data = [np.abs(data)/1e6 for data in all_data]

per_elec_dist_mean = []
per_elec_dist_max = []
per_elec_dist_min = []
for i in np.arange(len(fiber_diam)):
    #per_elec_dist_mean.append(all_data[i].nanmean(axis=(1,2)))
    #per_elec_dist_max.append(all_data[i].nanmax(axis=(1,2)))
    #per_elec_dist_min.append(all_data[i].nanmin(axis=(1,2)))
    per_elec_dist_mean.append(np.nanmean(all_data[i], axis=(1,2)))
    per_elec_dist_max.append(np.nanmax(all_data[i], axis=(1,2)))
    per_elec_dist_min.append(np.nanmin(all_data[i], axis=(1,2)))


f, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2, 3, sharey=True, sharex=True)
ax3.axis('off')
plotLocs = [ax1,ax2,ax4,ax5,ax6]
for j in np.arange(len(elec_dist)):
    curAx = plotLocs[j]
    for i in np.arange(len(fiber_diam)):
        curAx.plot(durations, per_elec_dist_mean[i][j], '-*k', color=edgecolors[i], label='%d um' % (fiber_diam[i]))
        curAx.fill_between(durations, per_elec_dist_min[i][j], per_elec_dist_max[i][j], alpha=0.5, edgecolor=edgecolors[i], facecolor=facecolors[i])
    curAx.axhline(y=0.5, color='r', linewidth=3)
    curAx.axhline(y=0.75, color='r', linewidth=3)
    curAx.set_xscale('log')
    curAx.set_xlabel('Duration (ms)')
    curAx.set_ylabel('Threshold current (mA)')
    curAx.set_title('elec_dist=%.2f mm' % (elec_dist[j]))
    if j==0:
        curAx.legend()
    curAx.set_ylim([0, 10])
    curAx.set_xlim([0.01, 2])
mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
plt.show(block=False)
plt.savefig('final.png',bbox_inches='tight')


"""
fiber_diam = 9
npzfile = np.load('threshCur_bipolarStim_%dum.npy' % (fiber_diam))

# For each diameter, want to average over elec_dist, node_offset, and dsep
# leaving only durations
npzfile = np.abs(npzfile)/1e6   # convert to mA
mean_val = npzfile.mean(axis=(0,1,2))
min_val = npzfile.min(axis=(0,1,2))
max_val = npzfile.max(axis=(0,1,2))

yerr = [mean_val-min_val, max_val-mean_val]
plt.figure()
plt.plot(durations, mean_val, '-*k', color='#CC4F1B')
plt.fill_between(durations, min_val, max_val, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
plt.xscale('log')
plt.xlabel('Duration (ms)')
plt.ylabel('Threshold current (mA)')
plt.title('%d um' % (fiber_diam))
plt.show(block=False)

#----
fiber_diam = 7
npzfile2 = np.load('threshCur_bipolarStim_%dum.npy' % (fiber_diam))
npzfile2 = np.abs(npzfile2)/1e6   # convert to mA
mean_val2 = npzfile2.mean(axis=(0,1,2))
min_val2 = npzfile2.min(axis=(0,1,2))
max_val2 = npzfile2.max(axis=(0,1,2))
yerr2 = [mean_val2-min_val2, max_val2-mean_val2]
plt.plot(durations, mean_val2, '-*k', color='#1B2ACC')
plt.fill_between(durations, min_val2, max_val, alpha=0.5, edgecolor='#1B2ACC', facecolor='#089FFF')
plt.show(block=False)

fiber_diam = 5
npzfile3 = np.load('threshCur_bipolarStim_%dum.npy' % (fiber_diam))
npzfile3 = np.abs(npzfile3)/1e6   # convert to mA
mean_val3 = npzfile3.mean(axis=(0,1,2))
min_val3 = npzfile3.min(axis=(0,1,2))
max_val3 = npzfile3.max(axis=(0,1,2))
yerr3 = [mean_val3-min_val3, max_val3-mean_val3]
plt.plot(durations, mean_val3, '-*k', color='#3F7F4C')
plt.fill_between(durations, min_val2, max_val, alpha=0.5, edgecolor='#3F7F4C', facecolor='#7EFF99')
plt.show(block=False)


#-------------------- Average for each elec_dist, durations
mean_at_dist = npzfile.mean(axis=(1,2))
plt.figure()
for i in np.arange(len(mean_at_dist)):
    plt.plot(durations, mean_at_dist[i], '-*r')
plt.xscale('log')
plt.xlabel('Duration (ms)')
plt.ylabel('Threshold current (mA)')

mean_at_dist2 = npzfile2.mean(axis=(1,2))
for i in np.arange(len(mean_at_dist2)):
    plt.plot(durations, mean_at_dist2[i], '-*b')

plt.show(block=False)
"""
