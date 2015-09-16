from matplotlib.pyplot import *
from colloids.phase import *
from os.path import isfile

qR = 0.1
q = qR2q(qR)

fig = figure('experimental phase diagram')
clf()
for phase, color in [('gel', 'r'), ('fluid', 'b'), ('transient', 'g'), ('clusters', 'y')]:
    phi, cp = np.loadtxt('phase_diag_%s.csv'%phase, unpack=True, skiprows=1, usecols=[1,2])
    sample = np.loadtxt('phase_diag_%s.csv'%phase, usecols=[0], skiprows=1, dtype='str')
    scatter(phi, cp, c=color, label=phase)
    for s, p, c in zip(np.atleast_1d(sample), np.atleast_1d(phi), np.atleast_1d(cp)):
        gca().annotate(s, (p,c))

#load or compute theoretical phase diagram
if isfile('phasediag.txt'):
    GL = np.loadtxt('phasediag.txt')
else:
    GL = Liu().all_GL(q, maxpiv=3500)
    np.savetxt('phasediag.txt', GL)
#draw binodal and spinodals
for f,l in zip(GL.T[1:], ['--']*2+['-']*2):
    cp = piv2y(GL[:,0], q) * alpha(f, qR)
    plot(100*f2vf(f), cp, 'k'+l)
#draw tie lines
for piv, fG, fL in GL[[50,60,70,80,85,88]+range(90,98),:3]:
    cpG = piv2y(piv, q) * alpha(fG, qR)
    cpL = piv2y(piv, q) * alpha(fL, qR)
    plot(100*f2vf(np.array([fG,fL])), [cpG, cpL], '-', color=(0.5,0.5,0.5))

#draw equicomposition lines
for x in np.arange(1,10)/10.:
    phiG, phiL = 100*f2vf(GL[:,1:3]).T
    phi = x*phiG+(1-x)*phiL
    cp = piv2y(GL[:,0], q) * alpha(vf2f(phi/100.), qR)
    plot(phi, cp, '-', color=(0.5,0.5,0.5))


ylabel(r'$c_p$ (g/L)')
xlabel(r'$\phi$ (%)')
ylim(0,3.5)
xlim(0,50)
draw()
savefig('phasediag.pdf')

fig = figure('theoretical phase diagram')
clf()

def xp2th(cp, phi, qR=0.1):
    """Converts experimental polymer concentration in osmotic pressure"""
    return y2piv(cp/alpha(vf2f(phi), qR2q(qR)), qR)

for phase, color in [('gel', 'r'), ('fluid', 'b'), ('transient', 'g'), ('clusters', 'y')]:
    phi, cp = np.loadtxt('phase_diag_%s.csv'%phase, unpack=True, skiprows=1, usecols=[1,2])
    sample = np.loadtxt('phase_diag_%s.csv'%phase, usecols=[0], skiprows=1, dtype='str')
    scatter(phi, xp2th(cp, phi*1e-2, qR), c=color, label=phase)
    for s, p, c in zip(np.atleast_1d(sample), np.atleast_1d(phi), np.atleast_1d(xp2th(cp, phi*1e-2, 0.1))):
        gca().annotate(s, (p,c))

#draw binodal and spinodals
for f,l in zip(GL.T[1:], ['--']*2+['-']*2):
    plot(100*f2vf(f), GL[:,0], 'k'+l)

#draw equicomposition lines
for x in np.arange(1,10)/10.:
    phiG, phiL = 100*f2vf(GL[:,1:3]).T
    plot(x*phiG+(1-x)*phiL, GL[:,0], '-', color=(0.5,0.5,0.5))

#for sb, l in zip('sb', ['-','--']):
 #   for ph in 'gl':
  #      phi, cp = np.loadtxt('gasliquid_%s%s.phd'%(sb,ph), unpack=True)
   #     plot(phi*100, xp2th(cp, phi, 0.1), l+'k')
ylabel(r'$\Pi v$')
xlabel(r'$\phi$ (%)')
ylim(0,3500)
xlim(0,74)
draw()
savefig('phasediagth.pdf')

