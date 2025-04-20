import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from numpy.polynomial.legendre import legval
from scipy.special import legendre

def get_gll_points(P):
    """Compute the GLL points (including endpoints)"""
    if P == 1:
        return np.array([-1.0, 1.0])
    Lp = legendre(P - 1)
    Lp_der = np.polyder(Lp)
    interior = np.roots(Lp_der)
    gll_points = np.concatenate(([-1.0], np.sort(interior), [1.0]))

    return gll_points
# Recursive relation for Jacobi-polynomials
def jacobi(n, x, a, b):
    a1 = 2 * n * (n + a + b) * (2*n + a + b - 2)
    a2 = (2*n - 1 + a + b) * (a**2 - b**2)
    a3 = (2*n - 2 + a + b) * (2*n - 1 + a + b) * (2*n + a + b)
    a4 = 2 * (n + a - 1) * (n + b - 1) * (2*n + a + b)
    
    if (n>1):
        return ( (a2 + a3*x) * jacobi(n-1, x, a, b) - a4 * jacobi(n-2, x, a, b) ) / a1
    elif (n==1):
        return 1 / 2 * ( a - b + (a + b + 2) * x)
    elif (n==0):
        return 1

def getmodified(xi, P, alpha, beta):

    nQ = len(xi)
    B = np.zeros((nQ, P+1))
    print(np.shape(B))

    for j in range(0,P+1):
        if j==0:
            for i in range(nQ):
                B[i,j] = (1-xi[i]) / 2
        elif j==P:
            for i in range(nQ):
                B[i,j] = (1+xi[i]) / 2
        else:
            for i in range(nQ):
                B[i,j] = (1-xi[i]**2) / 4  * jacobi(j-1, xi[i], alpha, beta)

        print(j)
    
    return B

def getlagrange(xi, P):

    """Compute the basis matrix B of shape (len(xi), P+1)"""
    xj = get_gll_points(P)  # GLL interpolation nodes
    nQ = len(xi)
    B = np.ones((nQ, P + 1))

    for j in range(P + 1):
        for i in range(nQ):
            # Compute L_j(xi[i])
            num = 1.0
            den = 1.0
            for k in range(P + 1):
                if k != j:
                    num *= (xi[i] - xj[k])
                    den *= (xj[j] - xj[k])
            B[i, j] = num / den
    return B

P = 4
Q = 4

fig = plt.figure()

gs = GridSpec(Q+2, P+2, figure=fig)

axlist = [[None] * (P+2) for _ in range(Q+2)]

for col in range(P+2):
    for row in range(Q+2):
        axlist[row][col] = fig.add_subplot(gs[row, col])


alpha = 1
beta = 1

nx = 11

xi = np.linspace(-1,1,nx)

B = getmodified(xi, P, alpha, beta)

for i in range(P+1):
    axlist[0][i+1].plot(xi, B[:,i]/np.max(abs(B[:,i])), label=r'$\psi_{'+str(i)+'}$'+r'$(\xi)$', linewidth=1.5)
    axlist[0][i+1].set_xlabel(r'$\xi_1$')
    axlist[0][i+1].set_title(fr'$\psi_{i}(\xi_1)$')
    axlist[0][i+1].grid()
    axlist[0][i+1].set_xlim([-1,1])
    axlist[0][i+1].set_ylim([-1,1])

for i in range(Q+1):
    axlist[i+1][0].plot(xi, B[:,i]/np.max(abs(B[:,i])), color='C1', label=r'$\psi_{'+str(i)+'}$'+r'$(\xi)$', linewidth=1.5)
    axlist[i+1][0].set_xlabel(r'$\xi_2$')
    axlist[i+1][0].set_ylabel(fr'$\psi_{i}(\xi_2)$', labelpad=10)
    axlist[i+1][0].grid()
    axlist[i+1][0].set_xlim([-1,1])
    axlist[i+1][0].set_ylim([-1,1])


for i in range(P+1):

    for j in range(Q+1):
        
        # b_tensor = np.outer(B[:,i]/np.max(abs(B[:,i])), B[:,j]/np.max(abs(B[:,j])))
        b_tensor = np.outer(B[:,j]/np.max(abs(B[:,j])), B[:,i]/np.max(abs(B[:,i])))
        im = axlist[i+1][j+1].imshow(np.rot90(b_tensor), extent=[-1, 1, -1, 1], interpolation='quadric', cmap='RdYlBu', vmin=-1, vmax=1)
        axlist[i+1][j+1].set_xlabel(r'$\xi_1$')
        axlist[i+1][j+1].set_ylabel(r'$\xi_2$')
        axlist[i+1][j+1].set_xlim([-1,1])
        axlist[i+1][j+1].set_ylim([-1,1])

axlist[0][0].set_axis_off()
cbar_ax = fig.add_axes([0.95, 0.12, 0.01, 0.6])  # [left, bottom, width, height
fig.colorbar(im, cax=cbar_ax, label=fr'$\psi_p(\xi_1)\psi(\xi_2)$')

fig.set_size_inches(9,9)
fig.subplots_adjust(wspace=1, hspace=1)
fig.savefig('modifiedBasis.pdf', dpi=300, format='pdf', bbox_inches='tight')
