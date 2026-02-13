from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
from matplotlib.patches import Circle
 
# first create create a function that evaluates the right-hand-side of the
# state-space equations for a given state vector

# x is the current state vector 
# t is current simulation time


def fibonacci(n):
    if n <= 1:
        return n
    return fibonacci(n - 1) + fibonacci(n - 2)

def fib_seq(n):

    fib_arr = []

    for i in range(n):
        fib_arr.append(fibonacci(i))

    return fib_arr


def toymodel(q, t, Re, A):
    
    v, eta = q
    dqdt = [A[0,0] * v + A[0,1] * eta, A[1,0] * v + A[1,1] * eta]

    return dqdt

Re = 15

# Computing the eigenvectors
# A = np.array([[-1/Re, 0], [1, -2/Re]])
A = np.array([[-1/Re, -1],[0, - 2/Re]])

xLim = Re*4
yLim = Re*2
# next, define a grid of points at which we will show arrows
x0 = np.linspace(-xLim,xLim,15,endpoint=1)
x1 = np.linspace(-yLim,yLim,15,endpoint=1)

# create a grid
X0,X1=np.meshgrid(x0,x1)

# projections of the trajectory tangent vector 
dX0 = np.zeros(X0.shape)
dX1 = np.zeros(X1.shape)

shape1,shape2 = X1.shape

for indexShape1 in range(shape1):
    for indexShape2 in range(shape2):
        dxdtAtX=toymodel([X0[indexShape1,indexShape2],X1[indexShape1,indexShape2]], 0, Re, A)
        dX0[indexShape1,indexShape2]=dxdtAtX[0]
        dX1[indexShape1,indexShape2]=dxdtAtX[1]     

fig = plt.figure()
gs = GridSpec(1, 3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0:2])
ax2 = fig.add_subplot(gs[0, 2])

# Integrating in time
t = np.linspace(0, 30, 201)

for x0 in np.linspace(-xLim, xLim, 5, endpoint=1):
    for y0 in np.linspace(-yLim*1, yLim*1, 2, endpoint=1):
        sol = odeint(toymodel, [x0,y0], t, (Re,A))
        ax1.plot(sol[:,0], sol[:,1], 'r', lw=0.01)

for y0 in np.linspace(-yLim, yLim, 12, endpoint = 1):
    for x0 in np.linspace(-xLim, xLim, 2, endpoint=1):
        sol = odeint(toymodel, [x0,y0], t, (Re,A))
        ax1.plot(sol[:,0], sol[:,1], 'r', lw=0.01)

# Computing the eigenmodes
[D,V] = np.linalg.eig(A)

# Optimal initial conditions
[L,S,R] = np.linalg.svd(A)
L = L*xLim/4
t_arr = np.linspace(0, 100, 201)
sol = odeint(toymodel, L[:,1], t_arr, (Re,A))
ax1.plot(sol[:,0], sol[:,1], 'C0', lw=1.2)
ax2.plot(t_arr, np.sqrt(sol[:,0]**2 + sol[:,1]**2), lw=1.5)

sol_trunc = sol[fib_seq(11),:]
for i in range(len(sol_trunc)):
    # ax1.quiver(0, 0, sol[i, 0], sol[i, 1], scale_units='x')
    ax1.annotate("", xytext=(0, 0), xy=(sol_trunc[i, 0], sol_trunc[i, 1]), arrowprops=dict(arrowstyle="->", color='C0', lw=0.7))
# ax1.quiver(sol[:-1, 0], sol[:-1, 1], sol[1:, 0]-sol[:-1, 0], sol[1:, 1]-sol[:-1, 1], scale_units='xy', angles='xy', scale=1, color='C0')    

# Random initial condition
angle = 200
norm = np.sqrt(L[0,1]**2 + L[1,1]**2)
rad = np.pi / 180 * angle
rand = [np.cos(rad)*norm, np.sin(rad)*norm]
sol = odeint(toymodel, rand, t_arr, (Re,A))
ax1.plot(sol[:,0], sol[:,1], 'C2', lw=1.2)
ax2.plot(t_arr, np.sqrt(sol[:,0]**2 + sol[:,1]**2), 'C2', lw=1.5)

sol_trunc = sol[fib_seq(9)[::2],:]
for i in range(len(sol_trunc)):
    # ax1.quiver(0, 0, sol[i, 0], sol[i, 1], scale_units='x')
    ax1.annotate("", xytext=(0, 0), xy=(sol_trunc[i, 0], sol_trunc[i, 1]), arrowprops=dict(arrowstyle="->", color='C2', lw=0.7))
# ax1.quiver(sol[:-1, 0], sol[:-1, 1], sol[1:, 0]-sol[:-1, 0], sol[1:, 1]-sol[:-1, 1], scale_units='xy', angles='xy', scale=1, color='C0')    

init_energy = np.sqrt(sol[0,0]**2 + sol[0,1]**2)
ax2.plot([t_arr[0], t_arr[-1]], [init_energy, init_energy], 'k--')
ax2.text(t_arr[-1]-10, init_energy+0.5, fr'$||\mathbf{{q}}_0||_2 = {init_energy:.0f}$', ha = 'right', va='bottom')
ax2.set_xlabel(r'$t$')
ax2.set_ylabel(r'$||\mathbf{q}||_2$')
ax2.grid()
ax2.set_xlim([t_arr[0], t_arr[-1]])
ax2.set_ylim([0, 60])
ax2.set_title(r'$(b)$')


# Plotting magnitude of opt
center = (0.0, 0.0)
circle = Circle(center, norm, fill=0, linestyle='--', lw=0.5)
ax1.add_patch(circle)


# Plotting stable node
ax1.scatter(0, 0, 20, 'k',zorder=100)
V1 = V[:,0] / np.linalg.norm(V[:,0])
# ax1.annotate("", xytext=(0, 0), xy=(V[0, 0], V[1, 0]), arrowprops=dict(arrowstyle="<-", color='k'))
# ax1.annotate("", xytext=(0, 0), xy=(V[0, 1], V[1, 1]), arrowprops=dict(arrowstyle="<-", color='k'))
# ax1.quiver(0, 0, V1[0], V1[1], angles='xy', scale_units='xy', scale=1, color='k')
scalars = np.linspace(10, 60, 8)  # from -2 to 2
# Base points for the arrows
starts = np.outer(scalars, V1)
# Arrow directions (all same as eigenmode direction)
xx = np.full_like(scalars, V1[0])
yy = np.full_like(scalars, V1[1])
# Plot
ax1.quiver(starts[:, 0], starts[:, 1], -xx, -yy, angles='xy', scale_units='xy', scale=0.33, color='k', width = 0.003, zorder=99)
ax1.plot([0, starts[-1,0]], [0, starts[-1,1]], 'k', lw = 1, zorder=90)
ax1.quiver(-starts[:, 0], -starts[:, 1], xx, yy, angles='xy', scale_units='xy', scale=0.33, color='k', width = 0.003, zorder=99)
ax1.plot([0, -starts[-1,0]], [0, -starts[-1,1]], 'k', lw = 1, zorder=90)

ax1.text(starts[-2,0], starts[-2,1]-1, r'$\mathbf{x_1}$', ha='right', va='top', fontsize=10)

V2 = V[:,1] / np.linalg.norm(V[:,1])
# ax1.annotate("", xytext=(0, 0), xy=(V[0, 0], V[1, 0]), arrowprops=dict(arrowstyle="<-", color='k'))
# ax1.annotate("", xytext=(0, 0), xy=(V[0, 1], V[1, 1]), arrowprops=dict(arrowstyle="<-", color='k'))
# ax1.quiver(0, 0, V2[0], V2[1], angles='xy', scale_units='xy', scale=1, color='k')
# scalars = np.linspace(0, 80, 21)  # from -2 to 2
# Base points for the arrows
starts = np.outer(scalars, V2)
# Arrow directions (all same as eigenmode direction)
xx = np.full_like(scalars, V2[0])
yy = np.full_like(scalars, V2[1])
# Plot
ax1.quiver(starts[:, 0], starts[:, 1], -xx, -yy, angles='xy', scale_units='xy', scale=0.33, color='k', width = 0.003, zorder=99)
ax1.plot([0, starts[-1,0]], [0, starts[-1,1]], 'k', lw = 1, zorder=90)
ax1.quiver(-starts[:, 0], -starts[:, 1], xx, yy, angles='xy', scale_units='xy', scale=0.33, color='k', width = 0.003, zorder=99)
ax1.plot([0, -starts[-1,0]], [0, -starts[-1,1]], 'k', lw = 1, zorder=90)

ax1.text(starts[-5,0], starts[-5,1]+1, r'$\mathbf{x_2}$', ha='right', va='bottom', fontsize=10)
# labelling eigenmodes



# ax1.quiver(X0, X1, dX0, dX1, scale=10, scale_units='xy')
# ax1.arrow(D[0]*V[0,0], D[0]*V[0,1], -D[0]*V[0,0], -D[0]*V[0,1])
# ax1.quiver(0,0, V[0,0], V[0,1], scale=1, scale_units='xy')
# ax1.plot([0, 10*V[0,0]], [0, 10*V[0,1]], zorder=99)
# ax1.plot([0, 10*V[1,0]], [0, 10*V[1,1]], zorder=99)
ax1.plot
ax1.grid()
ax1.set_xlim(-xLim*1., xLim*1.)
ax1.set_ylim(-yLim*1., yLim*1.)
ax1.set_xlabel(r'$v$')
ax1.set_ylabel(r'$\eta$')
ax1.set_title(r'$(a)$')

fig.set_size_inches(9, 3)
fig.subplots_adjust(wspace=0.3)
fig.savefig('PhasePotrait.pdf', dpi=300, bbox_inches='tight')
