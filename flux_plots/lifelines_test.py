from lifelines import KaplanMeierFitter
from lifelines.datasets import load_waltons
import matplotlib.pyplot as plt

df = load_waltons()

T = df['T']
E = df['E']

kmf = KaplanMeierFitter()
kmf.fit(T, E)

fig = plt.figure()
ax = kmf.plot()
plt.savefig('plots/KM_lifelines_test.png')
