import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from astropy.modeling import models, fitting 


#first - recreate example from online
#second - recreate example using 2d sample data
#third - use actual data
#fourth - clean up
#fifth - make work for multiple sources - loop thru

#creating random data
np.random.seed(42)
g1 = models.Gaussian2D(70, 100, 80, 40, 20, theta=0) #first gaussian to fit
g2 = models.Gaussian2D(80, 100, 100, 30, 20, theta=2) #second gaussian
x,y = np.mgrid[:128, :128]
z = g1(x,y) + g2(x,y)
z += np.random.normal(0, 5, z.shape) #adding noise to fake data

gg_init = models.Gaussian2D(76, 103, 76, 50, 15) + models.Gaussian2D(73, 102, 97, 27, 24) #input guess for fit
fitter = fitting.LevMarLSQFitter()
gg_fit = fitter(gg_init, x, y, z)

# Plot the data with the best-fit model
plt.figure(figsize=(8,3))
plt.subplot(1,3,1)
plt.imshow(z, origin='lower', interpolation='nearest', vmin=1e-2, vmax=120)
plt.title('data')

plt.subplot(1,3,2)
plt.imshow(gg_fit(x,y,z), origin='lower', interpolation='nearest', vmin=1e-2, vmax=120)
plt.title('model')

plt.subplot(1,3,3)
plt.imshow(z - gg_fit(x,y), origin='lower', interpolation='nearest', vmin=1e-2, vmax=120)
plt.title('residual')

plt.savefig('test_fitting.png')

#first - load in data, mask out to only get nearby region - take from gaussfit catalog
