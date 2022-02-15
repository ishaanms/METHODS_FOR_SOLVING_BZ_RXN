import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ

# Dimensionless parameters
ef=0.00012
e=0.036
f=1
q=0.00024

#packing 3 initial conditions with state of x,y,z into y0
initial_values = [1,1,1]
#Differential equation to be computed wrapped in a function
def Calculate (t,Y):
  x = Y[0]
  y = Y[1]
  z = Y[2]
  dxdt = (q*y - x*y+ x- x*x)*(1/e)
  dydt = (-1*q*y-x*y+f*z)*(1/ef)
  dzdt = x-z
  return [dxdt,dydt,dzdt]

#Computing the values
solution=integ.solve_ivp(Calculate,[0, 40],initial_values,method='RK45', atol=1e-4, rtol=1e-6)
Y=solution.y
timee=solution.t

#Scaling the result obtained for better visualization
xvalnew,yvalnew,zvalnew=[],[],[]
tempp=0
while tempp<len(Y[0]):
    xvalnew.append(np.log10(Y[0][tempp]))
    yvalnew.append(np.log10(Y[1][tempp]))
    zvalnew.append(np.log10(Y[2][tempp]))
    tempp+=1
#Plotting the graph
plt.plot(timee,xvalnew)        #Plotting HBrO2
plt.plot(timee,yvalnew)        #Plotting Br-
plt.plot(timee,zvalnew)        #Plotting Ce4+
plt.grid()
plt.show()
