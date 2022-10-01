import pyEGB
import matplotlib.pyplot as plt
import numpy as np
import math

size = 50 #number of models to be computed
minPressure = 0.1 #(x10^35 dyn/cm^2)
maxPressure = 30.0 #(x10^35 dyn/cm^2)
pressures = np.logspace(math.log10(minDensity),math.log10(maxDensity),size)
eosName = "ppwff1.cold.rns1.1.txt" #see EoS/ for other equations of state available
coupling = 20.0 # km^2
accuracy = 1e-5 # iterative method tolerance
maxIter = 500 # maximum method iterations
relaxation = 0.2 # relaxation factor
mass,radius, relaxationReturn, convergence = [np.zeros(size) for i in range(5)]
for centralPressure,i in zip(pressures,range(size)):
    mass[i],radius[i],relaxationReturn[i],convergence[i] = pyEGB.get_MR(eosName,coupling,centralPressure,accuracy,maxIter,relaxation)
    print(centralPressure,mass[i],radius[i])
    #print(centralPressure,mass[i],radius[i],relaxationReturn[i],convergence[i])

# plt.plot(radius,mass)
# plt.show()
