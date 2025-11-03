import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('FutureWork/sdc.csv', on_bad_lines='skip', header=0, names=['au', 'eps', 'pr', 'G', 'dom', 'sdc', 'rems'])

df = df.iloc[:,0:-1]

# Sampling Pr ~ 1
df = df.query('pr <= 1.2')
# Dropping white spaces
df['dom'] = df['dom'].apply(lambda x: x.replace(' ', ''))
df['sdc'] = df['sdc'].apply(lambda x: x.replace(' ', ''))

fig, ax = plt.subplots()
for index, item in df.iterrows():
    if item['dom'] == 'c':
        if item['sdc'] == 'Yes':
            ax.scatter(item['G'], item['eps'], color='C0', marker='o', edgecolor='k')
            ax.text(item['G']-2, item['eps']+0.015, item['au'], fontsize=8)
        elif item['sdc'] == 'No':
            ax.scatter(item['G'], item['eps'], color='C3', marker='o', edgecolor='k')
            ax.text(item['G']-2, item['eps']+0.015, item['au'], fontsize=8)
        else:
            ax.scatter(item['G'], item['eps'], color='C1', marker='o', edgecolor='k')
            ax.text(item['G']-2, item['eps']+0.015, item['au'], fontsize=8)
    elif item['dom'] == 'sq':
        if item['sdc'] == 'Yes':
            ax.scatter(item['G'], item['eps'], color='C0', marker='s', edgecolor='k')
            ax.text(item['G']-2, item['eps']+0.015, item['au'], fontsize=8)
        elif item['sdc'] == 'No':
            ax.scatter(item['G'], item['eps'], color='C3', marker='s', edgecolor='k')
            ax.text(item['G']-2, item['eps']+0.015, item['au'], fontsize=8)
        else:
            ax.scatter(item['G'], item['eps'], color='C1', marker='s', edgecolor='k')
            ax.text(item['G']-2, item['eps']+0.015, item['au'], fontsize=8)

# Our experiment
gammaRange = [6.28, 12.57, 25.13]
epsRange = [0.7, 0.7, 0.7]
color = ['C3', 'C1', 'C0']
for i in range(len(gammaRange)):
    G = gammaRange[i]
    EPS = epsRange[i]
    ax.scatter(G, EPS, color=color[i], marker='^', label=r'$\Gamma = $' + str(gammaRange[i]))

# Initial sample size
gammaRange = np.arange(10,60+10,10)
epsRange = np.arange(0.4,0.8+0.1,0.1)

for i in range(len(gammaRange)):
    G = gammaRange[i] * np.ones(len(epsRange))
    EPS = epsRange 
    ax.scatter(G, EPS, color='k', marker='x')

# Intended simulation
gammaRange = [20,30,40,50,60]
epsRange = [0.8, 0.7, 0.6, 0.5, 0.4]

for i in range(len(gammaRange)):
    G = gammaRange[i]
    EPS = epsRange[i]
    ax.scatter(G, EPS, color='C2', marker='x')

dfboden = pd.read_csv('Figures/bodenschatzfig31a.csv', names=['G', 'eps'])
ax.plot(dfboden['G'], dfboden['eps'], '-o')
ax.legend()

ax.set_ylabel(r'$\varepsilon$')
ax.set_xlabel(r'$\Gamma$')
ax.grid()

fig.set_size_inches([6,6])
fig.savefig('Figures/sdcompiled.pdf', format='pdf', dpi=300)
