from matplotlib.pyplot import *

fig = figure('phase diagram')
clf()
for phase, color in [('gel', 'r'), ('fluid', 'b'), ('transient', 'g')]:
    phi, cp = np.loadtxt('phase_diag_%s.csv'%phase, unpack=True, skiprows=1, usecols=[1,2])
    sample = np.loadtxt('phase_diag_%s.csv'%phase, usecols=[0], skiprows=1, dtype='str')
    scatter(phi, cp, c=color, label=phase)
    for s, p, c in zip(np.atleast_1d(sample), np.atleast_1d(phi), np.atleast_1d(cp)):
        gca().annotate(s, (p,c))

phi, cp = np.loadtxt('gasliquid_sg.phd', unpack=True)
plot(phi*100, cp, 'k')
ylabel(r'$c_p$ (g/L)')
xlabel(r'$\phi$ (%)')
ylim(0,3.5)
xlim(0,50)
draw()
savefig('phasediag.pdf')