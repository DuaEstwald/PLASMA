import numpy as np
import matplotlib.pyplot as plt


# =========================================================
#                       APARTADO 1
# =========================================================

# Utilizamos la ecuacion 

def dphi(w):
    dz = 1600*1e3 #m
    gamma = 5./3.
    g = 274 #m/s2
    H = 8.42*1e3 #m
    c = np.sqrt(gamma*g*H) #m/s
    wac = 3.7 #mH
    return (dz/c)*(w**2-wac**2)**(0.5)

w = np.linspace(0.1,10,100)

plt.figure()
plt.plot(w,dphi(w)*1e-1*np.pi/180,'k-.')
plt.ylabel(r'$\Delta\phi$ (rad)')
plt.xlabel('Freq (mHz)')
plt.savefig('fig1.png')
# ==========================================================
#                     APARTADO 2
# ==========================================================

import scipy.special as sp
z = np.linspace (0,50,10000)

# Utilizando la ecuacion 101 de A2:

#Supondremos 

# establecemos x,t = 0
t = 0
w=1.
H=1.
B=1.
rho_c=1.
mu0 = 4.*np.pi
A = 1.
V=B/(np.sqrt(mu0*rho_c*np.exp(-z/H)))
k=1.
x = 0.
h2 = (A*w*B/V)*sp.jv(1, 2*w*H/V)*np.exp(1j*(k*x+w*t))

plt.figure()

plt.plot (z,h2.real,'k-') 
plt.ylabel(r'$B_{1y}$')
plt.xlabel('h [km]')
plt.savefig('fig2.png')





# ========================================================
#                       APARTADO 3
# ========================================================

e = -1.6022e-19 #C
me = 9.10938e-31 #kg 
e0 = 8.8542e-12 #C2/Nm2
mH = 1.6735575e-27 #kg
kB = 1.3806485e-23 #J/K
sig_e = 1e-19 #m2
sig_i = 5e-19 #m2
cte_e = 3.7e-6
cte_i = 6e-8

mui = mH*mH/(mH+mH)
mue = me*mH/(me+mH)

ind,h,m,t5000,T,v,nH,ne,Pfotal,Pgas,sigma =np.loadtxt('tabla12.txt',unpack=True)

ne[np.where(ne>nH)] = nH[np.where(ne>nH)]
ni = ne
nn = nH-ni


L = 23.4-1.15*np.log10(ne)+3.45*np.log10(T)
def wpe(n):
    n_SI = n*((1e2)**3)
    return (n_SI*e**2/(me*e0))**0.5 #1/s


def wc(m,z,q):
    B = 100.*np.exp(-z/600.)*1e-4
    return np.abs(q)*B/m

def vn(n,T,sig,mu):
    n_SI = n*((1e2)**3)
    return n_SI*sig*(8*kB*T/(np.pi*mu))**0.5

def vie(n,T):
    n_SI = n*((1e2)**3)
    return (((e**4)*n_SI*L)/(3*(e0*me)**2))*(me/(2*np.pi*kB*T))**(3./2.)


def vii(n,T):
    n_SI = n*((1e2)**3)
    return (((e**4)*n_SI*L)/(3*np.sqrt(2)*(e0*mH)**2))*(mH/(2*np.pi*kB*T))**(3./2.)
#    return cte*n_SI*L/T**(3./2.)

plt.ion()

plt.figure()
plt.semilogy(h,wpe(ne),'k',label='Langmuir')
plt.semilogy(h,wc(mH,h,-e),'k:',label='Cyclotron (i)')
plt.semilogy(h,wc(me,h,e),'k--',label='Cyclotron (e)')
plt.semilogy(h,vn(nH,T,sig_i,mui),'r:',label='(in)') 
plt.semilogy(h,vn(nH,T,sig_e,mue),'r--',label='(en)') 
plt.semilogy(h,vie(ni,T),'b--',label='ei')
plt.semilogy(h,vii(ne,T),'b:',label='ii')
plt.xlabel('HEIGHT [KM]')
plt.ylabel(r'FREQUENCY [S$^{-1}$]')
plt.ylim(1e2,1e13)
plt.legend(loc=3,ncol=2)
plt.savefig('fig31.png')



def rd(T,n):
    n_SI = n*((1e2)**3)
    return (kB*T*e0/(n_SI*e**2))**(0.5)

def rc(z,m,T,q):
    B = 100.*np.exp(-z/600.)*1e-4
    return (2*kB*T*m)**(0.5)/(np.abs(q)*B)


plt.figure()
plt.semilogy(h,rd(T,ne),'k',label='Debye')
plt.semilogy(h,rc(h,mH,T,-e),'k:',label='Larmor i')
plt.semilogy(h,rc(h,me,T,e),'k--',label='Larmor e')
plt.xlabel('HEIGHT [KM]')
plt.ylabel('SPATIAL SCALE [m]')
plt.legend(loc=2)
plt.savefig('fig32.png')


