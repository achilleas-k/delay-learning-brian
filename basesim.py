#!/usr/bin/env python
from __future__ import print_function
from brian import *
from brian.compartments import *
from brian.library.ionic_currents import *
from brian.library.synapses import *


defaultclock.dt = dt = 0.1*ms
duration = 300*ms

length = 1000*um
nseg = 50
seg_length = length/nseg
ci = 1*uF/cm**2
gij = 0.02 * usiemens  # WAT IS THIS - PUT IT SOMEWHERE
diam = 2*um
radius = diam/2
area = pi*diam*seg_length
Ci = ci*area
e_leak = -70*mV
rMi = 12*kohm*cm**2
rL = 200*ohm*cm
Ri = rMi/area
Qi = seg_length/(pi*radius**2)
Rij = rL*Qi
g_pas = 0.00004*siemens/cm**2*area
e_pas = -70*mV

#print("Time constant =", Cm / gl)
#print("Space constant =", .5 * (diam / (gl * Ri)) ** .5)
print("Time constant = %s" % (Ri*Ci))

print("Setting up cable segments ...")
somaseg = nseg-1
equations = []
for i in range(nseg):
    equations += Equations(
        "dV/dt = (Iin+coupling+leak+pas)/Ci : volt",
        V="V_%i" % i,
        Iin="Iin_%i" % i,
        coupling="coupling_%i" % i,
        Ci=Ci,
        leak="leak_%i" % i,
        pas="pas_%i" % i,
    )
    leak_channel_eqn = "leak = (e_leak-V)/Ri : amp"
    equations += Equations(leak_channel_eqn,
                           leak="leak_%i" % i,
                           V="V_%i" % i,
                           Ri=Ri,
                           e_leak=e_leak)
    passive_channel_eqn = "pas = g_pas*(e_pas-V) : amp"
    equations += Equations(passive_channel_eqn,
                           pas="pas_%i" % i,
                           V="V_%i" % i,
                           g_pas=g_pas,
                           e_pas=e_pas)
    equations += Equations("Iin : amp", Iin="Iin_%i" % i)
    coupling_eqn = []
    if (i > 0):
        coupling_eqn.append("(V_pre-V_cur)/Rij")
    if (i < nseg-1):
        coupling_eqn.append("(V_next-V_cur)/Rij")
    if (i == nseg-1):
        coupling_eqn.append("(V_soma-V_cur)/Rij")
    coupling_eqn = "coupling = "+"+".join(coupling_eqn)+" : amp"
    equations += Equations(coupling_eqn,
                           coupling="coupling_%i" % i,
                           V_pre="V_%i" % (i-1),
                           V_next="V_%i" % (i+1),
                           V_cur="V_%i" % (i),
                           Rij=Rij)

    # TODO: Insert Na
equations += Equations(
    "dV/dt = (coupling+leak)/Ci : volt",
    V="V_soma",
    coupling="coupling_soma",
    Ci=Ci,
    leak="leak_soma"
)
equations += Equations("coupling = (V_pre-V_cur)/Rij : amp",
                       coupling="coupling_soma",
                       V_pre="V_%i" % (nseg-1),
                       V_cur="V_soma",
                       Rij=Rij)
equations += Equations("leak = (e_leak-V)/Ri : amp",
                       leak="leak_soma",
                       V="V_soma",
                       Ri=Ri,
                       e_leak=e_leak)

print("Setting up synapses ...")
synlocs = [int(nseg*rel) for rel in [0.1, 0.2, 0.3, 0.9]]
print("Creating neuron group ...")
neuron = NeuronGroup(1, model=equations)
for i in range(nseg):
    setattr(neuron, "V_%i" % i, e_leak)
setattr(neuron, "V_soma", e_leak)
print("Creating input spikes ...")
inspikes = SpikeGeneratorGroup(2, [(0, 100*ms)])
inconn6 = Connection(inspikes[0], neuron, state="V_%i" % synlocs[1])
inconn6.connect(inspikes, neuron, W=30*mV)

print("Creating monitors ...")
trace = {}
for i in synlocs+[nseg-1]:
    trace[i] = StateMonitor(neuron, "V_%i" % i, record=True)
trace["soma"] = StateMonitor(neuron, "V_soma", record=True)
spikemon = SpikeMonitor(neuron)
print("Running simulation for %s ..." % (duration))
run(duration)

for sli in synlocs:
    print("Segment %i peaked at %f ms with %f mV" % (
        sli, argmax(trace[sli][0])*dt*1000, (max(trace[sli][0])-float(e_leak))*1000))
print("Soma peaked at %f ms with %f mV" % (
    argmax(trace["soma"][0])*dt*1000, (max(trace["soma"][0])-float(e_leak))*1000))

figure()
subplot(211)
for sli in(synlocs):
    plot(trace[sli].times*1000, trace[sli][0]*1000, label="%i" % sli)
axis([0*second, duration*1000, -80, 0])
ylabel("voltage (mV)")
title("Dendrites")
legend(loc="best")
subplot(212)
plot(trace["soma"].times*1000, trace["soma"][0]*1000)
axis([0*second, duration*1000, -80, 0])
ylabel("voltage (mV)")
xlabel("time (ms)")
title("Soma")
#show()
