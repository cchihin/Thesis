import matplotlib.pyplot as plt
import numpy as np
import scipy.io

# Places text in the midpoint given by data
def annotate(data, label, ax, shiftx, shifty):

    mid_x = (max(data[:,0]) + min(data[:,0])) / 2
    mid_y = (max(data[:,1]) + min(data[:,1])) / 2
    print(mid_x)
    print(mid_y)
    ax.text(mid_x + shiftx , mid_y + shifty, label, va='center', ha='center')

osc = np.genfromtxt('Oscillatory.csv', delimiter=',')
sv = np.genfromtxt('SkewedVaricose.csv', delimiter=',')
eck = np.genfromtxt('Eckhuas.csv', delimiter=',')

# Sorting SV based on varepsilon
ind = np.argsort(sv[:,1])
sv = sv[ind,:]

# Sorting eck based on q
ind = np.argsort(eck[:,0])
eck = eck[ind,:]

# Reading neutral stability curve
neu_stab_x = scipy.io.loadmat('alpha_cr.mat')['x']
neu_stab_y = scipy.io.loadmat('eps_cr.mat')['y']

fig, ax = plt.subplots()

ax.plot(osc[:,0], osc[:,1])
ax.plot(sv[:,0], sv[:,1])
ax.plot(eck[:,0], eck[:,1])

ax.plot(neu_stab_x[0], neu_stab_y[0])

annotate(osc, 'Oscillatory', ax, -0.1, 0.5)
annotate(sv, 'Skewed-Varicose', ax, 0.05, 0)
annotate(eck, 'Eckhaus', ax, -0.8, -0.7)

ax.grid()
ax.set_xlabel(r'$q$')
ax.set_ylabel(r'$\varepsilon$')

fig.savefig('BusseBalloonPlapp.pdf')
