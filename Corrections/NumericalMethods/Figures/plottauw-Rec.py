import matplotlib.pyplot as plt
import numpy as np

def getretau(tau_w, nu, h):
    
    return np.sqrt(tau_w) * h / nu

def getdudy(u_c):

    return 2 * u_c

def Dean78_getturb(Re_c,U_c,rho=1):

    # Dean78's is based on U_b
    # We need to convert it to U_c
    C_f = 0.073 * (Re_c*4/3)**(-1/4)
    Ub = 2/3 * U_c
    tau_w_turb = C_f * 1/2 * rho * Ub**2

    return tau_w_turb

def Dean78_getlam(Re_c,U_c,rho=1):

    # Dean78's is based on U_b
    # We need to convert it to U_c
    Ub = 2/3 * U_c
    tau_w_lam = 12 / (Re_c*4/3) * 1/2 * rho * Ub**2

    return tau_w_lam

def readdata(fname, U_c, rho=1):

    Ub = 2./3 * U_c
    data = np.genfromtxt(f'{fname}', delimiter=',', skip_header=1)
    if data.ndim == 1:
        # x is based on Re_b -> converting to Re_c
        x = data[0] * 3/4
        # y is C_f -> converting to tau_w = Cf * 1/2 * rho * Ub**2
        y = data[1] * 1/2 * rho * Ub**2 
    else:
        # x is based on Re_b -> converting to Re_c
        x = data[:,0] * 3/4
        # y is C_f -> converting to tau_w = Cf * 1/2 * rho * Ub**2
        y = data[:,1] * 1/2 * rho * Ub**2 

    return x,y

U_c = 1
h = 1

# tau_w = np.arange(len(Re_c))
# tau_w = U_c*h/Re_c * getdudy(U_c)
# 
Re_ct = np.arange(2000, 5600, 5)
tauw_t = Dean78_getturb(Re_ct,U_c) * 2 / U_c / U_c

Re_cl = np.arange(400, 5600, 5)
tauw_l = Dean78_getlam(Re_cl,U_c) * 2 / U_c / U_c

fig, ax = plt.subplots()

# ax.plot(Re_c, tau_w, lw=3, label=r'Laminar flow')
ax.plot(Re_cl, tauw_l, lw=3, label=r'Laminar flow: $c_{f} = 4/Re_c$')
ax.plot(Re_ct, tauw_t, lw=3, label=r'Dean 1978: $c_{f} = 0.00302Re_c^{-1/4}$')

filelist = ['csvs/PatelHead1969.csv', 'csvs/Kim87.csv', 'csvs/IidaNagano1998.csv', 'csvs/Tsukahara2014.csv']
lbllist = ['Patel & Head 1969', 'Kim et. al. 1987', 'Iida & Nagano 1998', 'Tsukahara et. al. 2014']
msstyle = ['x', '*', '^', 'o']
collist = ['k', 'k', 'k', 'k']
for fname, ms, c, lbl in zip(filelist, msstyle, collist, lbllist):
    
    x, y = readdata(fname,U_c)
    # lbl = fname.split('.')[0].split('/')[1]
    ax.scatter(x, y*2, color=c, marker=ms, label=lbl, zorder=10)
    # ax[1].scatter(x, getretau(y, U_c*h/x, h), color=c, marker=ms, label=lbl, zorder=10)

    # if lbl == 'Kim87':
    #     y = getretau(y, U_c*h/x, h)
    #     # ax[1].text(x, y-5, f'({x:.0f}, {y:.1f})', va = 'top', ha = 'left')

ax.grid()
# ax.set_ylabel(r'$f = \tau_w = \nu \frac{dU}{dy} \frac{h}{U_c}$')
ax.set_ylabel(r'$c_f$')
ax.set_xlabel(r'$Re_c$')
# ax.set_title(r'$(a)$')
ax.set_ylim([0.0003*2, 0.003*1.9])
ax.set_xlim([500, 5500])
ax.legend(fontsize=9, loc=3)
# ax[1].grid()
# ax[1].set_ylabel(r'$Re_\tau = \frac{u_\tau h}{\nu}$')
# ax[1].set_xlabel(r'$Re_c$')
# ax[1].set_title(r'$(b)$')

# Plotting value of (Re_c, Re_tau) = (4200, 180)
# bf = 179.1**2 * 1/4200**2
# ax.plot([1100, 4200], [bf, bf], 'C3-.')
# ax.scatter(1100, bf, 100, color='C3', marker='X')
# ax.text(1200, bf*0.98, f'$(1100,{bf:.5f})$', ha='left', va='top')
# ax.plot([4200, 4200], [0, bf], 'C3-.')

fig.set_size_inches(8,4)
fig.savefig('cf-Rec.pdf', bbox_inches='tight', dpi=300)
