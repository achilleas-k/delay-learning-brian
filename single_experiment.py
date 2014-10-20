from __future__ import division
from brian import (
    NeuronGroup, SpikeGeneratorGroup, Synapses,
    EmpiricalThreshold, SpikeMonitor, StateMonitor,
    second, uF, msiemens, usiemens, nsiemens, mV, ms,
    run
)
import matplotlib.pyplot as plt

Cm = 1*uF # /cm**2
gL = 0.1*msiemens
EL = -65*mV
ENa = 55*mV
EK = -90*mV
gNa = 35*msiemens
gK = 9*msiemens
threshold = EmpiricalThreshold(threshold=15*mV, refractory=2*ms, state="V")

# Input parameters
taue = 5*ms
taui = 5*ms
EExc = 0*mV
EInh = -80*mV
WExc = 10*usiemens
WInh = 50*nsiemens

eqs='''
dV/dt=(-gNa*m**3*h*(V-ENa)\
    -gK*n**4*(V-EK)-gL*(V-EL)\
    -gExc*(V-EExc)\
    -gInh*(V-EInh)\
    +Iapp)/Cm : volt

m=alpham/(alpham+betam) : 1

alpham=-0.1/mV*(V+35*mV)/(exp(-0.1/mV*(V+35*mV))-1)/ms : Hz

betam=4*exp(-(V+60*mV)/(18*mV))/ms : Hz

dh/dt=5*(alphah*(1-h)-betah*h) : 1

alphah=0.07*exp(-(V+58*mV)/(20*mV))/ms : Hz

betah=1./(exp(-0.1/mV*(V+28*mV))+1)/ms : Hz

dn/dt=5*(alphan*(1-n)-betan*n) : 1

alphan=-0.01/mV*(V+34*mV)/(exp(-0.1/mV*(V+34*mV))-1)/ms : Hz

betan=0.125*exp(-(V+44*mV)/(80*mV))/ms : Hz

dgExc/dt = -gExc*(1./taue) : siemens

dgInh/dt = -gInh*(1./taui) : siemens

Iapp : amp

'''
neuron = NeuronGroup(1, eqs, threshold=threshold, method='RK')

neuron.V = -70*mV

# delays
B1, A1, A2, A3 = 5*ms, 20*ms, 30*ms, 40*ms
target_delay = A2-B1  # delay to be learned by neuron

spikes_A = [(0, 10*ms), (0, 115*ms), (0, 300*ms), (0, 450*ms)]
spikes_B = [(1, 10*ms), (1, 130*ms), (1, 335*ms), (1, 475*ms)]
inputs = SpikeGeneratorGroup(2, spikes_A+spikes_B)
synapse_A = Synapses(inputs[0], neuron,
                     model="w : siemens", pre="gExc_post += w")
synapse_A[:,:] = 3
synapse_A.w = WExc
synapse_A.delay[0] = A1
synapse_A.delay[1] = A2
synapse_A.delay[2] = A3

synapse_B = Synapses(inputs[1], neuron,
                     model="w : siemens", pre="gExc_post += w")
synapse_B[:,:] = 1
synapse_B.w = WExc
synapse_B.delay[0] = B1



vtrace = StateMonitor(neuron, 'V', record=True)
mtrace = StateMonitor(neuron, 'm', record=True)
htrace = StateMonitor(neuron, 'h', record=True)
ntrace = StateMonitor(neuron, 'n', record=True)
exc_trace = StateMonitor(neuron, 'gExc', record=True)

spike_monitor = SpikeMonitor(neuron)

run(1.0*second)

print("Number of spikes fired: %i" % (spike_monitor.nspikes))
if spike_monitor.nspikes <= 10:
    print("Spike times follow:")
    print(', '.join([str(sp) for sp in spike_monitor[0]]))


plt.subplot(211)
plt.plot(vtrace.times, vtrace[0], label="V(t)")
plt.legend()
#subplot(312)
#plot(mtrace.times, mtrace[0], label="m(t)")
#plot(htrace.times, htrace[0], label="h(t)")
#plot(ntrace.times, ntrace[0], label="n(t)")
#legend()
plt.subplot(212)
plt.plot(exc_trace.times, exc_trace[0], label="gExc")
plt.legend()
plt.show()

