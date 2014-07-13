import numpy as np
import matplotlib.pyplot as plt

## Load in data
TS = np.load("../Data/Data_Ts.npy")
plt.plot(TS)
plt.show()


