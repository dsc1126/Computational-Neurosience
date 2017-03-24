import matplotlib.pyplot as plt
import math
import numpy as np
from random import randint


# Q1
# Input parameters
# Time period from 0 to 1 second
Tau_m = 10 * 10**-3		# 10 ms
EL = Vr = -70 * 10**-3	# -70 mV	
Vt = -40 * 10**-3		# -40 mV
Rm = 10 * 10**6			# 10 mega ohm
Ie = 3.1 * 10**-9		# 3.1 nA
dt = 1 * 10**-3			# 0.001 s time interval

# Define differential equation
def integrateFire( V, t ) :
    return ( EL - V + Rm*Ie ) / Tau_m

# Simulate neuron
t0=0
t1=1
dt=0.001

V0=Vr

ts=[t0]
Vs=[V0]

while ts[-1]<t1:
    
    V=Vs[-1]+integrateFire(Vs[-1],ts[-1])*dt
    
    ts.append(ts[-1]+dt)
    Vs.append(V)
    if Vs[-1] > Vt :	# If exceed threshold, make the value as Vr
        Vs[-1] = Vr

plt.figure(1)
plt.plot(ts,Vs)
plt.xlabel("time (in Seconds)")
plt.ylabel("Voltage (in Volts)")
plt.axis([ 0, 1, -0.08, -0.03 ])
plt.title("Q1 Simulate an Integrate and Fire Model")
plt.savefig("Q1.png", dpi=500)
plt.show()


# Q2
# Calculate minimum current Ie to produce an action potential
Iemin = (Vt - EL) / Rm 	# Iemin = 3.0 nA
print("minimum current Ie required is",Iemin)	#display result in cmd window


# Q3
# Assume Ie is 0.1 nA lower than Ie minmum value required from Q2
Ie = Iemin - 0.1 * 10 ** -9		# Ie = 2.9 nA

# Input parameters
# Time period from 0 to 1 second
Tau_m = 10 * 10**-3		# 10 ms
EL = Vr = -70 * 10**-3	# -70 mV	
Vt = -40 * 10**-3		# -40 mV
Rm = 10 * 10**6			# 10 mega ohm
dt = 1 * 10**-3			# 0.001 s time interval

# Define differential equation
def integrateFire( V, t ) :
    return ( EL - V + Rm*Ie ) / Tau_m

# Simulate neuron
t0=0
t1=1
dt=0.001

V0=Vr

ts=[t0]
Vs=[V0]

while ts[-1]<t1:
    
    V=Vs[-1]+integrateFire(Vs[-1],ts[-1])*dt
    
    ts.append(ts[-1]+dt)
    Vs.append(V)
    if Vs[-1] > Vt :	# If exceed threshold, make the value as Vr
        Vs[-1] = Vr

plt.figure(3)
plt.plot(ts,Vs)
plt.xlabel("time (in Seconds)")
plt.ylabel("Voltage (in Volts)")
plt.axis([ 0, 1, -0.08, -0.03 ])
plt.title("Q3 Simulate an Integrate and Fire Model (with Ie = 2.9 nA)")
plt.savefig("Q3.png", dpi=500)
plt.show()


# Q4
# find number of spikes in 1 s in relation with the function of currents ranging from 2 to 5 nA
Ies = np.arange(2,5.1,step=0.1) * 10**-9		# from 2 to 5 nA
# Input parameters
# Time period from 0 to 1 second
Tau_m = 10 * 10**-3		# 10 ms
EL = Vr = -70 * 10**-3	# -70 mV	
Vt = -40 * 10**-3		# -40 mV
Rm = 10 * 10**6			# 10 mega ohm
dt = 1 * 10**-3			# 0.001 s time interval

# Define differential equation
def integrateFire( V, t ) :
    return ( EL - V + Rm*Ie ) / Tau_m

# Simulate neuron
t0=0
t1=1
dt=0.001

V0=Vr

ts=[t0]
Vs=[V0]
spike_num = []

for Ie in Ies:
	ts=[t0]
	Vs=[V0]
	# Ie = c
	spikes = 0
	while ts[-1]<t1:
		V=Vs[-1]+integrateFire(Vs[-1],ts[-1])*dt
		ts.append(ts[-1]+dt)
		Vs.append(V)
		if Vs[-1] > Vt :	# If exceed threshold, make the value as Vr
			Vs[-1] = Vr
			spikes += 1
	spike_num.append(spikes)

plt.figure(4)
plt.plot(Ies,spike_num)
plt.xlabel("Current Ie (in Amperes)")
plt.ylabel("Number of Spikes")
plt.axis([ 1.5*10**-9, 5.5*10**-9, -10, 120 ])
plt.title("Q4 Spiking Rate in relation with Current Ie")
plt.savefig("Q4.png", dpi=500)
plt.show()


# Q5
# Input parameters
Tau_m = 20 * 10 ** -3		# 20 ms
EL = -70 * 10 ** -3 		# -70 mV
Vr = -80 * 10 ** -3 		# -80 mV
Vt = -54 * 10 ** -3 		# -54 mV
RmIe = 18 * 10 ** -3 		# 18 mV
RmGs = 0.15            		# 0.15
P = 0.5            			# 0.5
Tau_s = 10 * 10 ** -3 		# 10 ms	
dt = 1 * 10**-3				# 0.001 s time interval

# Define differential equation
def integrateFire2( V, t, tf ) :
    Ist = RmGs*(P*math.exp(-(t-tf)/Tau_s))*(Es-V)	#single exponential synapse model
    return (EL-V+RmIe+Ist)/Tau_m

# Q5a
Es = 0	#excitatory synapse
VaInitial = randint(Vr * 10 **3, Vt * 10 **3) * 10 ** -3 # randomly between Vr and Vt
VbInitial = randint(Vr * 10 **3, Vt * 10 **3) * 10 ** -3 # randomly between Vr and Vt

t0=0
t1=1
dt=0.001
ts=[t0]

Va = []
Va.append(VaInitial)
ta = 0

Vb = []
Vb.append(VbInitial)
tb = 0

while ts[-1]<t1:
    
	V1 = Va[-1]+dt*integrateFire2(Va[-1],ts[-1]-dt,tb) 
	V2 = Vb[-1]+dt*integrateFire2(Vb[-1],ts[-1]-dt,ta)
	Va.append(V1)
	Vb.append(V2)
	ts.append(ts[-1]+dt)
	if Va[-1] > Vt :
		Va[-1] = Vr
		ta = ts[-1]
	if Vb[-1] > Vt :
		Vb[-1] = Vr
		tb = ts[-1]	


# Print results
plt.figure(5)
p1, = plt.plot(ts,Va,'b')
p2, = plt.plot(ts,Vb,'r')
plt.xlabel("time (in Seconds)")
plt.ylabel("Voltage (in Volts)")
plt.title("Q5(a) Integrate and Fire Model \n of Two Connected Neurons (Excitatory Synapses)")
plt.legend([p1, p2], ["Neuron A", "Neuron B"])
plt.axis([ 0, 1, -0.085, -0.045 ])
plt.savefig("Q5a.png", dpi=500)
plt.show()

# Q5b
Es = Vr	#inhibitory synapse
VaInitial = randint(Vr * 10 **3, Vt * 10 **3) * 10 ** -3 # randomly between Vr and Vt
VbInitial = randint(Vr * 10 **3, Vt * 10 **3) * 10 ** -3 # randomly between Vr and Vt

t0=0
t1=1
dt=0.001
ts=[t0]

Va = []
Va.append(VaInitial)
ta = 0

Vb = []
Vb.append(VbInitial)
tb = 0

while ts[-1]<t1:
    
	V1 = Va[-1]+dt*integrateFire2(Va[-1],ts[-1]-dt,tb) 
	V2 = Vb[-1]+dt*integrateFire2(Vb[-1],ts[-1]-dt,ta)
	Va.append(V1)
	Vb.append(V2)
	ts.append(ts[-1]+dt)
	if Va[-1] > Vt :
		Va[-1] = Vr
		ta = ts[-1]
	if Vb[-1] > Vt :
		Vb[-1] = Vr
		tb = ts[-1]	


# Print results
plt.figure(6)
p1, = plt.plot(ts,Va,'b')
p2, = plt.plot(ts,Vb,'r')
plt.xlabel("time (in Seconds)")
plt.ylabel("Voltage (in Volts)")
plt.title("Q5(b) Integrate and Fire Model \n of Two Connected Neurons (Inhibitory Synapses)")
plt.legend([p1, p2], ["Neuron A", "Neuron B"])
plt.axis([ 0, 1, -0.085, -0.045 ])
plt.savefig("Q5b.png", dpi=500)
plt.show()


# Q6
# Input parameters
# Time period from 0 to 1 second
Tau_m = 10 * 10**-3		# 10 ms
EL = Vr = -70 * 10**-3	# -70 mV	
Vt = -40 * 10**-3		# -40 mV
Rm = 10 * 10**6			# 10 mega ohm
Ie = 3.1 * 10**-9		# 3.1 nA
dt = 1 * 10**-3			# 0.001 s time interval
Ek = -80 * 10**-3		# -80 mV
Delta_Gm = 0.0005* 10**-6# 0.005 mega ohm ^-1 conductance
#Delta_Rm = 1 / Delta_Gm # find the resistance change rate
decay_Tau = 200 * 10**-3# 200 ms	

# Define differential equation
def integrateFire3( V, t ) :
    return ( EL - V + Rm*Ie + Rm*(Gk+Gk_delta)*(Ek - V)) / Tau_m

# Simulate neuron
Gk=0
Gk_delta=0
t0=0
t1=1
dt=0.001
V0=Vr

ts=[t0]
Vs=[V0]

while ts[-1]<t1:

    V=Vs[-1]+integrateFire3(Vs[-1],ts[-1])*dt
    Gk += -(dt*Gk)/decay_Tau
    ts.append(ts[-1]+dt)
    Vs.append(V)
    if Vs[-1] > Vt :	# If exceed threshold, make the value as Vr
        Vs[-1] = Vr
        Gk = 0
        Gk_delta += Delta_Gm  
    

plt.figure(7)
plt.plot(ts,Vs)
plt.xlabel("time (in Seconds)")
plt.ylabel("Voltage (in Volts)")
plt.axis([ 0, 1, -0.08, -0.03 ])
plt.title("Q6 Simulate an Integrate and Fire Model with Potassium Current")
plt.savefig("Q6.png", dpi=500)
plt.show()


