#!/usr/bin/env python
'''
Dendrite with 100 compartments
'''
from __future__ import print_function
from brian import *
from brian.compartments import *
from brian.library.ionic_currents import *
from brian.library.synapses import *

defaultclock.dt = dt = 0.1*ms
duration = 100*ms

length = 1000 * um
nseg = 50
dx = length / nseg
Cm = 1 * uF / cm ** 2
gl = 0.02 * msiemens / cm ** 2
diam = 12 * um
area = pi * diam * dx
El = 0 * mV
Ri = 330 * ohm * cm
ra = Ri * 4 / (pi * diam ** 2)
tau_mem = 10*ms

#print("Time constant =", Cm / gl)
#print("Space constant =", .5 * (diam / (gl * Ri)) ** .5)

print("Setting up cable segments ...")
segments = {}
somaseg = nseg-1
for i in range(nseg):
    segments[i] = MembraneEquation(Cm * area)
    segments[i] += leak_current(gl * area, El)

#segments[0] += Current("I : amp")
print("Setting up synapses ...")
synlocs = [int(nseg*rel) for rel in [0.1, 0.2, 0.3, 0.9]]
for sl in synlocs:
    segments[i] += exp_synapse(
        input='Iin', tau=10*ms, unit=amp, output='__membrane_Im')
cable = Compartments(segments)
for i in range(1, nseg):
    cable.connect(i-1, i, ra*dx)

print("Creating neuron group ...")
lastmem = "vm_"+str(somaseg)
neuron = NeuronGroup(1, model=cable)
#neuron.V_0=10*mV
#neuron.I_0 = .02 * nA
print("Creating input spikes ...")
inspikes = SpikeGeneratorGroup(2, [(0, 20*ms)])
inconn0 = Connection(inspikes[0], neuron, state="Iin_"+str(synlocs[1]))
inconn0.connect(inspikes, neuron, W=1*mV)
#inconn6 = Connection(inspikes[1], neuron, state="vm_6")
#inconn6.connect(inspikes, neuron, W=40*mV)

print("Creating monitors ...")
trace = []
for i in synlocs+[nseg-1]:
    trace.append(StateMonitor(neuron, 'vm_' + str(i), record=True))
spikemon = SpikeMonitor(neuron)
print("Running simulation for %s ..." % (duration))
run(duration)

for idx, sl in enumerate(synlocs):
    print("Segment %i peaks at %f ms with %f mV" % (
        sl, argmax(trace[idx][0])*dt*1000, max(trace[idx][0])*1000))
trace[-1].insert_spikes(spikemon, 40*mV)

for idx, sl in enumerate(synlocs):
    plot(trace[idx].times*1000, trace[idx][0], label="%i" % sl)
ylabel("voltage (V)")
legend(loc="best")
show()
