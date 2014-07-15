#!/usr/bin/env python
'''
Dendrite with 100 compartments
'''
from __future__ import print_function
from brian import *
from brian.compartments import *
from brian.library.ionic_currents import *
from brian.library.synapses import *

#@network_operation
#def current(clock):
#    if clock.t < 1*ms:
#        neuron.Iin_10 = 0.1*uA
#    else:
#        neuron.Iin_10 = 0*uA

defaultclock.dt = dt = 0.1*ms
duration = 100*ms

length = 1000 * um
nseg = 50
dx = length / nseg
ci = 1 * uF / cm**2
gij = 0.02 * usiemens
diam = 12 * um
radius = diam/2
area = pi * diam * dx
Ci = ci*area
El = 0 * mV
rMi = 10*Mohm*cm**2
rL = 330*ohm*cm
Ri = rMi/area
Qi = dx/(pi*radius**2)
Rij = rL*Qi

#print("Time constant =", Cm / gl)
#print("Space constant =", .5 * (diam / (gl * Ri)) ** .5)

print("Setting up cable segments ...")
somaseg = nseg-1
equations = []
for i in range(nseg):
    equations += Equations("dV/dt = (Iin-V/Ri+coupling)/Ci : volt",
                           V="V_%i" % i,
                           Iin="Iin_%i" % i,
                           coupling="coupling_%i" % i,
                           Ri=Ri,
                           Ci=Ci)
    equations += Equations("Iin : amp", Iin="Iin_%i" % i)
    if (i > 0) & (i < nseg-1):
        equations += Equations("coupling = ((V_pre-V_cur) + (V_post-V_cur))/Rij : amp",
                               coupling="coupling_%i" % i,
                               V_pre="V_%i" % (i-1),
                               V_post="V_%i" % (i+1),
                               V_cur="V_%i" % (i),
                               Rij=Rij)
    elif (i == 0):
        equations += Equations("coupling = (V_post-V_cur)/Rij : amp",
                               coupling="coupling_%i" % i,
                               V_post="V_%i" % (i+1),
                               V_cur="V_%i" % (i),
                               Rij=Rij)
    elif (i == nseg-1):
        equations += Equations("coupling = (V_pre-V_cur)/Rij : amp",
                               coupling="coupling_%i" % i,
                               V_pre="V_%i" % (i-1),
                               V_cur="V_%i" % (i),
                               Rij=Rij)

print("Setting up synapses ...")
synlocs = [int(nseg*rel) for rel in [0.1, 0.2, 0.3, 0.9]]
print("Creating neuron group ...")
neuron = NeuronGroup(1, model=equations)
neuron.V_10 = 0.1*mV
print("Creating input spikes ...")
inspikes = SpikeGeneratorGroup(2, [(0, 20*ms)])
#inconn0 = Connection(inspikes[0], neuron, state="Iin_"+str(synlocs[1]))
#inconn0.connect(inspikes, neuron, W=1*mV)
#inconn6 = Connection(inspikes[1], neuron, state="vm_6")
#inconn6.connect(inspikes, neuron, W=40*mV)

print("Creating monitors ...")
trace = []
coup = []
for i in synlocs+[nseg-1]:
    trace.append(StateMonitor(neuron, "V_%i" % i, record=True))
    coup.append(StateMonitor(neuron, "coupling_%i" % i, record=True))
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
