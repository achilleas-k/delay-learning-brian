from __future__ import print_function
from brian import *
from brian.compartments import *
from brian.library.ionic_currents import *
from brian.library.synapses import *

# TODO: Remove all channels and figure out if there's something wrong with
# calcium parameters or equations.
# TODO: Conductances per area! Add "/cm**2" to channel conductances and set
# channel current equation units to "amp/metre**2"
# TODO: Pay attention to differences between dendritic compartments and soma
# TODO: Major cleanup: Split code into functions or even files, simplify adding
# of channels (lots of repetition in channel code), make compartmentalisation
# easier to specify in the future without all this mess

def leak(seg, g, e):
    # leak channel
    leak_channel_eqn = "leak = g_leak*(e_leak-V) : amp"
    equations = Equations(leak_channel_eqn,
                          leak="leak_"+str(seg),
                          V="V_"+str(seg),
                          g_leak=g,
                          e_leak=e)
    return equations

def external(seg):
    # external current
    equations = Equations("Iin : amp", Iin="Iin_"+str(seg))
    return equations

def excitatory(seg, e, tau):
    # excitatory inputs (synapses)
    excitatory_equation = "excIn = g_exc*(e_exc-V) : amp"
    equations = Equations(excitatory_equation,
                           excIn="excIn_"+str(seg),
                           V="V_"+str(seg),
                           g_exc="g_exc_"+str(seg),
                           e_exc=e)
    exc_conductance_equation = "dg_exc/dt = -g_exc*(1./tau_exc) : siemens"
    equations += Equations(exc_conductance_equation,
                           dg_exc="dg_exc_"+str(seg),
                           g_exc="g_exc_"+str(seg),
                           tau_exc=tau)
    return equations

def sodium(seg, g, e):
    # sodium (Na)
    sodium_channel_eqn = "Na = g_Na*m**3*h*(e_Na-V) : amp"
    equations = Equations(sodium_channel_eqn,
                          Na="Na_"+str(seg),
                          m="m_"+str(seg),
                          h="h_"+str(seg),
                          V="V_"+str(seg),
                          g_Na=g,
                          e_Na=e)
    sodium_state_eqns = """
    m = alpha_m/(alpha_m+beta_m) : 1
    alpha_m = -0.1/mV*(V+35*mV)/(exp(-0.1/mV*(V+35*mV))-1) : 1
    beta_m = 4*exp(-(V+60*mV)/(18.0*mV)) : 1

    dh/dt = (5*(alpha_h*(1-h)-beta_h*h))/ms : 1
    alpha_h = 0.07*exp(-(V+58*mV)/(20*mV)) : 1
    beta_h = 1.0/(exp(-0.1/mV*(V+28*mV))+1) : 1
    """
    equations += Equations(sodium_state_eqns,
                           m="m_"+str(seg),
                           alpha_m="alpha_m_"+str(seg),
                           beta_m="beta_m_"+str(seg),
                           dh="dh_"+str(seg),
                           h="h_"+str(seg),
                           alpha_h="alpha_h_"+str(seg),
                           beta_h="beta_h_"+str(seg),
                           V="V_"+str(seg))
    return equations

def calcium_t(seg, g, e):
    # t-type calcium (CaT)
    calcium_t_channel_eqn = "CaT = g_CaT*r**3*s*(e_CaT-V) : amp"
    equations = Equations(calcium_t_channel_eqn,
                           CaT="Ca_"+str(seg),
                           r="r_"+str(seg),
                           s="s_"+str(seg),
                           e_CaT=e,
                           g_CaT=g,
                           V="V_"+str(seg))
    calcium_t_state_eqns = """
    dr/dt = (alpha_r*(1-r)-beta_r*r)/ms : 1
    alpha_r = 1.0/(1.7+exp(-(V+28.2*mV)/13.5)) : 1
    beta_r = exp(-(V+63.0*mV)/7.8)/(exp(-(V+28.8*mV)/13.1)+1.7) : 1

    ds/dt = (alpha_s*(1-s-d)-beta_s*s)/ms : 1
    alpha_s = exp(-(V+160.3*mV)/17.8) : 1
    beta_s  = ((0.25+exp((V+83.5*mV)/6.3))**0.5-0.5) * (exp(-(V+160.3*mV)/17.8)) : 1

    dd/dt = (beta_d*(1-s-d)-alpha_d*d)/ms : 1
    bd     = (0.25+exp((V+83.5*mV)/6.3))**0.5 : 1
    alpha_d = (1.0+exp((V+37.4*mV)/30.0))/(240.0*(0.5+bd)) : 1
    beta_d  = (bd-0.5)*alpha_d : 1
    """
    equations += Equations(calcium_t_state_eqns,
                           V="V_"+str(seg),
                           r="r_"+str(seg),
                           alpha_r="alpha_r_"+str(seg),
                           beta_r="beta_r_"+str(seg),
                           s="s_"+str(seg),
                           alpha_s="alpha_s_"+str(seg),
                           beta_s="beta_s_"+str(seg),
                           d="d_"+str(seg),
                           alpha_d="alpha_d_"+str(seg),
                           beta_d="beta_d_"+str(seg),
                           bd="bd_"+str(seg)
                           )
    return equations

def calcium_l(seg, g, e):
    # l-type calcium (CaL)
    # TODO: Adapt from calm95.mod (in neuron/delay_learning/channels/)
    pass


def potassium(seg, g, e):
    # potassium (K)
    potassium_channel_eqn = "K = g_K*n**4*(e_K-V) : amp"
    equations = Equations(potassium_channel_eqn,
                           K="K_"+str(seg),
                           n="n_"+str(seg),
                           e_K=e,
                           g_K=g,
                           V="V_"+str(seg))
    potassium_state_eqns = """
    dn/dt = 5*(alpha_n*(1-n)-beta_n*n)/ms : 1
    alpha_n = -0.01/mV*(V+34*mV)/(exp(-0.1/mV*(V+34*mV))-1) : 1
    beta_n = 0.125*exp(-(V+44*mV)/(80*mV)) : 1
    """
    equations += Equations(potassium_state_eqns,
                           dn="dn_"+str(seg),
                           n="n_"+str(seg),
                           alpha_n="alpha_n_"+str(seg),
                           beta_n="beta_n_"+str(seg),
                           V="V_"+str(seg))
    return equations

def coupling(seg, R, nseg):
    # coupling equations
    coupling_eqn = []
    if (seg > 0):
        coupling_eqn.append("(V_pre-V_cur)/Rij")
    if (seg < nseg-1):
        coupling_eqn.append("(V_next-V_cur)/Rij")
    if (seg == nseg-1):
        coupling_eqn.append("(V_soma-V_cur)/Rij")
    if coupling_eqn:
        coupling_eqn = "coupling = "+"+".join(coupling_eqn)+" : amp"
    else:
        coupling_eqn = "coupling : amp"
    equations = Equations(coupling_eqn,
                          coupling="coupling_"+str(seg),
                          V_pre="V_"+str(seg-1),
                          V_next="V_"+str(seg+1),
                          V_cur="V_"+str(seg),
                          Rij=R)
    return equations


defaultclock.dt = dt = 0.1*ms
duration = 300*ms

length = 1000*um  # total length
nseg = 1  # number of segments (dendrite)
seg_length = length/nseg  # segment length
ci = 1*uF/cm**2  # specific membrane capacitance
gij = 0.02 * usiemens  # WAT IS THIS - PUT IT SOMEWHERE
diam = 2*um  # diameter
radius = diam/2  # radius
seg_area = pi*diam*seg_length  # segment surface area
Ci = 6*uF  #ci*area  # membrane capacitance
rMi = 12*kohm*cm**2  # specific membrane resistance
rL = 200*ohm*cm  # intracellular/longitudinal resistivity
Ri = 12*kohm  #rMi/area  # membrane resistance
Qi = seg_length/(pi*radius**2)  # axial resistance factor
Rij = rL*Qi  # coupling resistance
Rij = Ri/10

## CHANNELS ##
g_leak = 0.3*msiemens  # leak conductance
e_leak = -70*mV  # membrane leak potential

g_Na = 35*msiemens  #/cm**2*seg_area  # sodium conductance
e_Na = 75*mV  # sodium reversal potential

g_K = 9*msiemens  # potassium conductance
e_K = -90*mV  # potassium reversal potential

g_CaT = 2*msiemens  # calcium t-type conductance
e_CaT = 130*mV  # calcium t-type reversal potential

g_CaL = 3*msiemens  # calcium l-type conductance


e_exc = 0*mV  # excitatory reversal potential
tau_exc = 15*ms  # excitatory conductance time constant
w_exc = 100*uS  # excitatory weight (conductance change)

#print("Time constant =", Cm / gl)
#print("Space constant =", .5 * (diam / (gl * Ri)) ** .5)
print("Time constant = %s" % (Ci/g_leak))

print("Setting up cable segments ...")
somaseg = nseg-1
equations = []
for i in range(nseg):
    equations += Equations(
        "dV/dt = (Iin+coupling+leak+Na+K+CaT+excIn)/Ci : volt",
        V="V_%i" % i,
        Iin="Iin_%i" % i,
        coupling="coupling_%i" % i,
        Ci=Ci,
        leak="leak_%i" % i,
        Na="Na_%i" % i,
        excIn="excIn_%i" % i,
        K="K_%i" % i,
        CaT="Ca_%i" % i,
    )
    equations += leak(i, g_leak, e_leak)
    equations += external(i)
    equations += excitatory(i, e_exc, tau_exc)
    equations += sodium(i, g_Na, e_Na)
    equations += calcium_t(i, g_CaT, e_CaT)
    equations += potassium(i, g_K, e_K)
    equations += coupling(i, Rij, nseg)

equations += Equations(
    "dV/dt = (coupling+Na+K+leak+CaT)/Ci : volt",
    V="V_soma",
    coupling="coupling_soma",
    Ci=Ci,
    leak="leak_soma",
    Na="Na_soma",
    K="K_soma",
    CaT="Ca_soma"
)
equations += leak("soma", g_leak, e_leak)
equations += sodium("soma", g_Na, e_Na)
equations += calcium_t("soma", g_CaT, e_CaT)
equations += potassium("soma", g_K, e_K)
equations += Equations("coupling = (V_pre-V_cur)/Rij : amp",
                       coupling="coupling_soma",
                       V_pre="V_%i" % (nseg-1),
                       V_cur="V_soma",
                       Rij=Rij)

print("Setting up synapses ...")
synlocs = [int(nseg*rel) for rel in [0.1, 0.2, 0.3, 0.6]]
print("Creating neuron group ...")
neuron = NeuronGroup(1, model=equations)
for i in range(nseg):
    setattr(neuron, "V_%i" % i, e_leak)
setattr(neuron, "V_soma", e_leak)
print("Creating input spikes ...")
inspikes = SpikeGeneratorGroup(2, [(0, 150*ms), (1, 38*ms)])
inconn_a2 = Connection(inspikes[0], neuron, state="g_exc_%i" % synlocs[1])
inconn_a2.connect(inspikes, neuron, W=w_exc)
#inconn_a3 = Connection(inspikes[0], neuron, state="g_exc_%i" % synlocs[2])
#inconn_a3.connect(inspikes, neuron, W=w_exc)
#inconn_b = Connection(inspikes[1], neuron, state="g_exc_%i" % synlocs[-1])
#inconn_b.connect(inspikes, neuron, W=w_exc)

print("Creating monitors ...")
trace = {}
states = []
for i in synlocs+[nseg-1]:
    trace[i] = StateMonitor(neuron, "V_%i" % i, record=True)
    states.extend(["Ca_%i" % i, "Na_%i" % i, "K_%i" % i])
channel_traces = MultiStateMonitor(neuron, states, record=True)
trace["soma"] = StateMonitor(neuron, "V_soma", record=True)
spikemon = SpikeMonitor(neuron)
print("Running simulation for %s ..." % (duration))
run(duration, report="stdout")

for sli in synlocs:
    print("Segment %i peaked at %f ms with %f mV" % (
        sli, argmax(trace[sli][0])*dt*1000, (max(trace[sli][0])-float(e_leak))*1000))
    print("\tVoltage sum: %f mV" % (sum(trace[sli][0]-min(trace[sli][0]))*1000))
print("Soma peaked at %f ms with %f mV" % (
    argmax(trace["soma"][0])*dt*1000, (max(trace["soma"][0])-float(e_leak))*1000))
print("\tVoltage sum: %f mV" % (sum(trace["soma"][0]-min(trace["soma"][0]))*1000))

ion()
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
