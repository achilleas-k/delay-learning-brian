#!/usr/bin/env python
from __future__ import print_function
from brian import *
from brian.compartments import *
from brian.library.ionic_currents import *
from brian.library.synapses import *


defaultclock.dt = dt = 0.1*ms
duration = 300*ms

length = 1000*um  # total length
nseg = 10  # number of segments (dendrite)
seg_length = length/nseg  # segment length
ci = 1*uF/cm**2  # specific membrane capacitance
gij = 0.02 * usiemens  # WAT IS THIS - PUT IT SOMEWHERE
diam = 2*um  # diameter
radius = diam/2  # radius
seg_area = pi*diam*seg_length  # segment surface area
Ci = 1*uF  #ci*area  # membrane capacitance
e_leak = -70*mV  # membrane leak potential
rMi = 12*kohm*cm**2  # specific membrane resistance
rL = 200*ohm*cm  # intracellular/longitudinal resistivity
Ri = 12*kohm  #rMi/area  # membrane resistance
Qi = seg_length/(pi*radius**2)  # axial resistance factor
Rij = rL*Qi  # coupling resistance
Rij = Ri
g_pas = 0.00004*siemens/cm**2*seg_area  # passive channel conductance
e_pas = -70*mV  # passive channel reversal potential
g_Na = 35*msiemens/cm**2*seg_area  # sodium conductance
e_Na = 55*mV  # sodium reversal potential
g_K = 9*msiemens  # potassium conductance
e_K = -90*mV  # potassium reversal potential
e_exc = 0*mV  # excitatory reversal potential
tau_exc = 15*ms  # excitatory conductance time constant
w_exc = 80*uS  # excitatory weight (conductance change)

#print("Time constant =", Cm / gl)
#print("Space constant =", .5 * (diam / (gl * Ri)) ** .5)
print("Time constant = %s" % (Ri*Ci))

print("Setting up cable segments ...")
somaseg = nseg-1
equations = []
for i in range(nseg):
    equations += Equations(
        "dV/dt = (Iin+coupling+leak+pas+Na+K+excIn)/Ci : volt",
        V="V_%i" % i,
        Iin="Iin_%i" % i,
        coupling="coupling_%i" % i,
        Ci=Ci,
        leak="leak_%i" % i,
        pas="pas_%i" % i,
        Na="Na_%i" % i,
        excIn="excIn_%i" % i,
        K="K_%i" % i,
    )
    # leak channel
    leak_channel_eqn = "leak = (e_leak-V)/Ri : amp"
    equations += Equations(leak_channel_eqn,
                           leak="leak_%i" % i,
                           V="V_%i" % i,
                           Ri=Ri,
                           e_leak=e_leak)
    # passive channel
    passive_channel_eqn = "pas = g_pas*(e_pas-V) : amp"
    equations += Equations(passive_channel_eqn,
                           pas="pas_%i" % i,
                           V="V_%i" % i,
                           g_pas=g_pas,
                           e_pas=e_pas)
    # external current
    equations += Equations("Iin : amp", Iin="Iin_%i" % i)
    # excitatory inputs (synapses)
    excitatory_equation = "excIn = g_exc*(e_exc-V) : amp"
    equations += Equations(excitatory_equation,
                           excIn="excIn_%i" % i,
                           V="V_%i" % i,
                           g_exc="g_exc_%i" % i,
                           e_exc=e_exc)
    exc_conductance_equation = "dg_exc/dt = -g_exc*(1./tau_exc) : siemens"
    equations += Equations(exc_conductance_equation,
                           dg_exc="dg_exc_%i" % i,
                           g_exc="g_exc_%i" % i,
                           tau_exc=tau_exc)
    # sodium (Na)
    sodium_channel_eqn = "Na = g_Na*m**3*h*(e_Na-V) : amp"
    equations += Equations(sodium_channel_eqn,
                           Na="Na_%i" % i,
                           m="m_%i" % i,
                           h="h_%i" % i,
                           V="V_%i" % i,
                           g_Na=g_Na,
                           e_Na=e_Na)
    equations += Equations("m = alpha_m/(alpha_m+beta_m) : 1",
                           m="m_%i" % i,
                           alpha_m="alpha_m_%i" % i,
                           beta_m="beta_m_%i" % i)
    equations += Equations("alpha_m = -0.1/mV*(V+35*mV)/(exp(-0.1/mV*(V+35*mV))-1)/ms : Hz",
                           alpha_m="alpha_m_%i" % i,
                           V="V_%i" % i)
    equations += Equations("beta_m = 4*exp(-(V+60*mV)/(18.0*mV))/ms : Hz",
                           beta_m="beta_m_%i" % i,
                           V="V_%i" % i)
    equations += Equations("dh/dt = 5*(alpha_h*(1-h)-beta_h*h) : 1",
                           dh="dh_%i" % i,
                           h="h_%i" % i,
                           alpha_h="alpha_h_%i" % i,
                           beta_h="beta_h_%i" % i)
    equations += Equations("alpha_h = 0.07*exp(-(V+58*mV)/(20*mV))/ms : Hz",
                           alpha_h="alpha_h_%i" % i,
                           V="V_%i" % i)
    equations += Equations("beta_h = 1.0/(exp(-0.1/mV*(V+28*mV))+1)/ms : Hz",
                           beta_h="beta_h_%i" % i,
                           V="V_%i" % i)
    # potassium (K)
    potassium_channel_eqn = "K = g_K*n**4*(e_K-V) : amp"
    equations += Equations(potassium_channel_eqn,
                           K="K_%i" % i,
                           n="n_%i" % i,
                           e_K=e_K,
                           g_K=g_K,
                           V="V_%i" % i)
    equations += Equations("dn/dt = 5*(alpha_n*(1-n)-beta_n*n) : 1",
                           dn="dn_%i" % i,
                           n="n_%i" % i,
                           alpha_n="alpha_n_%i" % i,
                           beta_n="beta_n_%i" % i)
    equations += Equations("alpha_n = -0.01/mV*(V+34*mV)/(exp(-0.1/mV*(V+34*mV))-1)/ms : Hz",
                           alpha_n="alpha_n_%i" % i,
                           V="V_%i" % i)
    equations += Equations("beta_n = 0.125*exp(-(V+44*mV)/(80*mV))/ms : Hz",
                           beta_n="beta_n_%i" % i,
                           V="V_%i" % i)


    # coupling equations
    coupling_eqn = []
    if (i > 0):
        coupling_eqn.append("(V_pre-V_cur)/Rij")
    if (i < nseg-1):
        coupling_eqn.append("(V_next-V_cur)/Rij")
    if (i == nseg-1):
        coupling_eqn.append("(V_soma-V_cur)/Rij")
    if coupling_eqn:
        coupling_eqn = "coupling = "+"+".join(coupling_eqn)+" : amp"
    else:
        coupling_eqn = "coupling : amp"
    equations += Equations(coupling_eqn,
                           coupling="coupling_%i" % i,
                           V_pre="V_%i" % (i-1),
                           V_next="V_%i" % (i+1),
                           V_cur="V_%i" % (i),
                           Rij=Rij)

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
synlocs = [int(nseg*rel) for rel in [0.1, 0.2, 0.3, 0.6]]
print("Creating neuron group ...")
neuron = NeuronGroup(1, model=equations)
for i in range(nseg):
    setattr(neuron, "V_%i" % i, e_leak)
setattr(neuron, "V_soma", e_leak)
print("Creating input spikes ...")
inspikes = SpikeGeneratorGroup(2, [(0, 100*ms)])
inconn6 = Connection(inspikes[0], neuron, state="g_exc_%i" % synlocs[-1])
inconn6.connect(inspikes, neuron, W=w_exc)

print("Creating monitors ...")
trace = {}
for i in synlocs+[nseg-1]:
    trace[i] = StateMonitor(neuron, "V_%i" % i, record=True)
trace["soma"] = StateMonitor(neuron, "V_soma", record=True)
spikemon = SpikeMonitor(neuron)
print("Running simulation for %s ..." % (duration))
run(duration, report="stdout")

for sli in synlocs:
    print("Segment %i peaked at %f ms with %f mV" % (
        sli, argmax(trace[sli][0])*dt*1000, (max(trace[sli][0])-float(e_leak))*1000))
print("Soma peaked at %f ms with %f mV" % (
    argmax(trace["soma"][0])*dt*1000, (max(trace["soma"][0])-float(e_leak))*1000))

figure()
subplot(211)
for sli in(synlocs):
    plot(trace[sli].times*1000, trace[sli][0]*1000, label="%i" % sli)
#axis([0*second, duration*1000, -80, 40])
ylabel("voltage (mV)")
title("Dendrites")
legend(loc="best")
subplot(212)
plot(trace["soma"].times*1000, trace["soma"][0]*1000)
#axis([0*second, duration*1000, -80, 0])
ylabel("voltage (mV)")
xlabel("time (ms)")
title("Soma")
#show()
