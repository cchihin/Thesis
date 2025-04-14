import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import scipy

def secondorder(f, x, ind):
    
    fpp = f[ind+1] - 2 * f[ind] + f[ind-1] / dxdx

alpha_cr = sio.loadmat('alpha_cr.mat')['x'][0]
eps_cr = sio.loadmat('eps_cr.mat')['y'][0]
ra_cr = eps_cr*1708.8 + 1708.8

ind_cr = np.argwhere(eps_cr == 0)[0,0]

# Approach 1 to computing gradients
dff = np.gradient(np.gradient(ra_cr, alpha_cr), alpha_cr)

# Approach 2 to computing gradients
spl = scipy.interpolate.splrep(alpha_cr,ra_cr,k=3) # no smoothing, 3rd order spline
ddy = scipy.interpolate.splev(alpha_cr,spl,der=2) # use those knots to get second derivative 
dddy = scipy.interpolate.splev(alpha_cr,spl,der=3) # use those knots to get second derivative 

g2 = ddy[ind_cr] * 1/2
g3 = dddy[ind_cr] * 1/3
print(g2)

# Import critical curve from plapp
crit = np.genfromtxt('pri0.71crit.csv', delimiter=',')

eckplapp = np.genfromtxt('pr071Plapp.csv', delimiter=',')

eckdecker = np.genfromtxt('EckhausDecker94.csv', delimiter=',')

# Pr = 1.07
# L = 0.17374 + 0.07792/Pr + 0.12927/Pr/Pr
# print(L)
# coef = (1 + L)/L

y = 3 * g2 * (alpha_cr - 3.117)**2  + 5 * g3 * (alpha_cr - 3.117)**3 + 1708.8

eqn1415 = np.genfromtxt('eqn1415.csv', delimiter=',')
eqn18 = np.genfromtxt('eqn18.csv', delimiter=',')

eck = np.genfromtxt('Eckhuas.csv', delimiter=',')

fig, ax = plt.subplots()

ax.plot(alpha_cr, ra_cr)
ax.scatter([1.87,4.375],[1.7*1708.8, 1.7*1708.8])

ax.plot(eck[:,0], (eck[:,1]*1708.8)+1708.8)

ax.plot(alpha_cr.flatten(), y.flatten())
ax.plot(crit[:,0], 10**crit[:,1])
ax.plot(eckplapp[:,0], 10**eckplapp[:,1])
# ax.plot(eqn1415[:,0], eqn1415[:,1])
# ax.plot(eqn18[:,0], eqn18[:,1])
ax.plot(eckdecker[:,0], (eckdecker[:,1]*1708.8)+1708.8)
ax.grid()
ax.set_ylim([0, 4*1708])
ax.set_xlim([0.5, 6])

fig.savefig('EckhausInstability.pdf', bbox_inches='tight', dpi=100)

