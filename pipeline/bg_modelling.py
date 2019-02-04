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
g1 = models.Gaussian2D(5, 1, 2, 0.4, 0.2, theta=0) #first gaussian to fit
g2 = models.Gaussian2D(3, 0.5, 0.5, 0.3, 0.2, theta=2) #second gaussian
x,y = np.mgrid[:128, :128]
z = g1(x,y) + g2(x,y)
z += np.random.normal(0, 0.1, z.shape) #adding noise to fake data
    
gg_init = models.Gaussian2D(4, 0.5, 2.3, 0.3, 0.3) + models.Gaussian2D(2, 0.5, 0.8, 0.4, 0.3) #input guess for fit
fitter = fitting.SLSQPLSQFitter()
gg_fit = fitter(gg_init, x, y, z)

# Plot the data with the best-fit model
plt.figure(figsize=(8,5))
plt.plot(x, y, 'ko')
plt.plot(x, gg_fit(x))
plt.xlabel('Position')
plt.ylabel('Flux')
plt.savefig('test_fitting.png')

#first - load in data, mask out to only get nearby region - take from gaussfit catalog
