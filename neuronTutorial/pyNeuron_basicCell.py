"""
Following http://neuron.yale.edu/neuron/static/docs/neuronpython/ballandstick1.html
Builds two-compartment cell - a soma and dendrite.
"""
from neuron import h, gui


"""
Sections - cell sections
"""
soma = h.Section(name='soma')
dend = h.Section(name='dend')

"""
Topology - connectivity between sections
"""
dend.connect(soma(1))

"""
Geometry - 3D location of the sections
"""
# Surface area of cylinder is 2*pi*r*h (sealed ends are implicit).
# Here we make a square cylinder in that the diameter
# is equal to the height, so diam = h. ==> Area = 4*pi*r^2
# We want a soma of 500 microns squared:
# r^2 = 500/(4*pi) ==> r = 6.2078, diam = 12.6157
soma.L = soma.diam = 12.6157
dend.L = 200    # microns
dend.diam = 1   # microns
print "Surface area of soma=", h.area(0.5, sec=soma)

shape_window = h.PlotShape()
shape_window.exec_menu('Show Diam')

"""
Biophysics - ionic channels and membrane properties of the sections
"""
for sec in h.allsec():
    sec.Ra = 100
    sec.cm = 1

# Active HH current in soma
soma.insert('hh')
soma.gnabar_hh = 0.12
soma.gkbar_hh = 0.036
soma.gl_hh = 0.0003
soma.el_hh = -54.3

# Passive current in dendrite
dend.insert('pas')
dend.g_pas = 0.001
dend.e_pas = -65

"""
Synapses - list of synapses onto the cell
"""

"""
Instrumentation
"""
stim = h.IClamp(dend(1))
stim.delay = 5  # starts 5ms after simulation starts
stim.dur = 1
stim.amp = 0.1

"""
Record vec, run and plot
"""
v_vec = h.Vector()
t_vec = h.Vector()
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)
simdur = 25.0

h.tstop = simdur
h.run()

# Initial plot of soma
from matplotlib import pyplot
pyplot.figure(figsize=(8,4))
pyplot.plot(t_vec, v_vec)
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.title('Soma voltage, stim amp=%s'%(stim.amp))
pyplot.show(block=False)

# Vary current amplitude in loop
import numpy
pyplot.figure(figsize=(8,4))
step = 0.075
num_steps = 4
for i in numpy.linspace(step, step*num_steps, num_steps):
    stim.amp = i
    h.tstop = simdur
    h.run()
    pyplot.plot(t_vec, v_vec, color='black')

pyplot.xlabel('time (ms')
pyplot.ylabel('mV')
pyplot.title('Soma voltage, varying amp')
pyplot.show(block=False)

# Also record voltage at dendrite
dend_v_vec = h.Vector()
dend_v_vec.record(dend(0.5)._ref_v)
pyplot.figure(figsize=(8,4))
for i in numpy.linspace(step, step*num_steps, num_steps):
    stim.amp = i
    h.tstop = simdur
    h.run()
    soma_plot = pyplot.plot(t_vec, v_vec, color='black')
    dend_plot = pyplot.plot(t_vec, dend_v_vec, color='red')
pyplot.legend(soma_plot + dend_plot, ['soma','dend'])
pyplot.xlabel('time (ms')
pyplot.ylabel('mV')
pyplot.title('Soma and dendrite voltage, varying amp')
pyplot.show(block=False)

# Check the effects of nseg in dendrite on the signal through dendrite
pyplot.figure(figsize=(8,4))
for i in numpy.linspace(step, step*num_steps, num_steps):
    stim.amp = i
    h.run()
    soma_plot = pyplot.plot(t_vec, v_vec, color='black')
    dend_plot = pyplot.plot(t_vec, dend_v_vec, color='red')
# change nseg then do it again
dend.nseg = 101
for i in numpy.linspace(step, step*num_steps, num_steps):
    stim.amp = i
    h.run()
    soma_hires= pyplot.plot(t_vec, v_vec, color='blue')
    dend_hires= pyplot.plot(t_vec, dend_v_vec, color='green')

pyplot.legend(soma_plot + dend_plot + soma_hires + dend_hires, ['soma','dend', 'soma hi-res', 'dend hi-res'])
pyplot.xlabel('time (ms')
pyplot.ylabel('mV')
pyplot.title('Soma and dendrite vs hi-res voltage, varying amp')
pyplot.show(block=False)

# Use even higher nseg, see how the error drops off
pyplot.figure(figsize=(8,4))
ref_v = []
ref_dend_v = []

# Run through to get high res values
dend.nseg = 101
for i in numpy.linspace(step, step*num_steps, num_steps):
    stim.amp = i
    h.run()
    soma_hires= pyplot.plot(t_vec, v_vec, color='blue')
    dend_hires= pyplot.plot(t_vec, dend_v_vec, color='green')
    # copy into reference vectors for comparison
    ref_v_vec = numpy.zeros_like(v_vec)
    v_vec.to_python(ref_v_vec)
    ref_v.append(ref_v_vec)

    ref_dend_v_vec = numpy.zeros_like(dend_v_vec)
    dend_v_vec.to_python(ref_dend_v_vec)
    ref_dend_v.append(ref_dend_v_vec)

# Run through low res
dend.nseg = 1   # start here, use odd values
err = 0
idx = 0
for i in numpy.arange(step, step*(num_steps+0.9), step):
    stim.amp = i
    h.run()
    soma_lowres = pyplot.plot(t_vec, v_vec, color='black')
    dend_lowres = pyplot.plot(t_vec, dend_v_vec, color='red')
    err += numpy.mean(numpy.abs(numpy.subtract(ref_v[idx], v_vec)))
    err += numpy.mean(numpy.abs(numpy.subtract(ref_dend_v[idx], dend_v_vec)))
    idx += 1

err /= idx
err /= 2 # Since we have a soma and dend vec
print "Average error =", err



