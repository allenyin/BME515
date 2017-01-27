from neuron import h, gui   # h for neuron interaction, gui to run simulation
from matplotlib import pyplot   # For plotting results
import pickle

soma = h.Section(name='soma')
soma.insert('pas')                  # passive channel in soma
asyn = h.AlphaSynapse(soma(0.5))    # insert an alpha synapse to soma segment
asyn.onset = 20
asyn.gmax = 1

# recording variables
v_vec = h.Vector()
t_vec = h.Vector()
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)
h.tstop = 40.0
h.run()

# plotting
pyplot.figure(figsize=(8,4))
pyplot.plot(t_vec, v_vec)
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show(block=False)

# save the vectors from simulation
with open('t_vec.p', 'w') as t_vec_file:
    pickle.dump(t_vec.to_python(), t_vec_file)
with open('v_vec.p', 'w') as v_vec_file:
    pickle.dump(v_vec.to_python(), v_vec_file)

def loadAndPlot():
    with open('t_vec.p') as t_vec_file:
        py_t_vec = pickle.load(t_vec_file)
    #t_vec_restore = h.Vector(py_t_vec)
    with open('v_vec.p') as v_vec_file:
        py_v_vec = pickle.load(v_vec_file)
    #v_vec_restore = h.Vector(py_v_vec)
    pyplot.figure(figsize=(8,4))
    pyplot.plot(py_t_vec, py_v_vec)
    pyplot.xlabel('time (mS)')
    pyplot.ylabel('mV')
    pyplot.title('Imposter!')
    pyplot.show(block=False)
