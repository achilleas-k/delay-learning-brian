#!/usr/bin/env python
'''
Dendrite with 100 compartments
'''
from __future__ import print_function
from brian import *
from brian.compartments import *
from brian.library.ionic_currents import *

length = 1 * mm
nseg = 10
dx = length / nseg
Cm = 1 * uF / cm ** 2
gl = 0.02 * msiemens / cm ** 2
diam = 1 * um
area = pi * diam * dx
El = 0 * mV
Ri = 100 * ohm * cm
ra = Ri * 4 / (pi * diam ** 2)
tau_mem = 10*ms

#print("Time constant =", Cm / gl)
#print("Space constant =", .5 * (diam / (gl * Ri)) ** .5)

segments = {}
segsies = {}
for i in range(nseg):
    segments[i] = MembraneEquation(Cm * area) + leak_current(gl * area, El)
#segments[0] += Current("I : amp")

cable = Compartments(segments)
for i in range(1, nseg):
    cable.connect(i-1, i, ra*dx)

lastmem = "vm_"+str(nseg-1)
neuron = NeuronGroup(1, model=cable, threshold=lastmem+">=20*mV",
                     reset=lastmem+"=0*mV")
#neuron.V_0=10*mV
#neuron.I_0 = .02 * nA
inspikes = SpikeGeneratorGroup(2, [(0, 20*ms), (1, 60*ms)])
inconn0 = Connection(inspikes[0], neuron, state="vm_0")
inconn0.connect(inspikes, neuron, W=40*mV)
#inconn6 = Connection(inspikes[1], neuron, state="vm_6")
#inconn6.connect(inspikes, neuron, W=40*mV)

trace = []
for i in range(nseg):
    trace.append(StateMonitor(neuron, 'vm_' + str(i), record=True))
spikemon = SpikeMonitor(neuron)
run(200 * ms)
trace[-1].insert_spikes(spikemon, 40*mV)
subplot(211)
for i in [6]:
    plot(trace[i].times, trace[i][0], label="%i" % i)
ylabel("voltage (V)")
legend(loc="best")
subplot(212)
plot(trace[-1].times, trace[-1][0])
xlabel("time (s)")
ylabel("voltage (V)")
show()
