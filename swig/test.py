import pyEGB
import matplotlib.pyplot as plt

mass,radius,relaxationReturn,convergence = pyEGB.get_MR("ppsly4.cold.rns1.1.txt", 30.0, 25.0, 1e-5,500,0.2)

print(mass,radius,relaxationReturn,convergence)