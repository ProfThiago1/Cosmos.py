from scipy.integrate import odeint
import numpy as np
import math
import matplotlib.pylab as plt
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
import pandas as pd

#Constantes
pc = 3.0857e16# m                                                                                               
kpc = 1e3 * pc
Mpc = 1e3 * kpc
Msun = 1.98847e30
G = 6.674184e-11  # Nm^2/kg^2
Pi = math.pi

#Função de Hubble H0(h)
def Ho(h):
    return 100*h*1000/Mpc #Pelo o que entendi, 1000 -> km

def rhoC(h):
    return 3*Ho(h)**2/(8*Pi*G)

def rho_de(a,h,om,w):
    return rhoC(h)*(1-om)*a**(-3*(1+w))

def rho_m(a,h,om):
    return rhoC(h)*om*a**(-3)

#Função parâmetro de densidade para M
def OmegaM(a,om,w):
    return 1/(1 +((1-om)/om)*a**(-3*w))

#Função parâmetro de densidade para E
def OmegaE(a,om,w):
    return 1 - OmegaM(a,om,w)

#Função densidade Energia
def densE(a,om,w):
    return OmegaE(a,om,w)*((3*Ho(0.7)**2)/(8*math.pi*G))*(om*a**(-3) + (1-om))

#função mtar sem pertubação em DE
def mtar(m,h,om,w):  
    return (-3*m*Msun/(4*Pi*rho_de(1,h,om,w)*(1+3*w)))**(1/3) /Mpc

def mtar_popolo(m,h,om,w):  
    return (-3*m*Msun/( 4*Pi*rho_de(1,h,om,w)*(1+3*w) - 4*Pi*rho_m(1,h,om)*alfa(m)*delta_c_m(m) ))**(1/3) /Mpc

#Função mtar com perturbação em DE
def mtar_de(m,h,om,w):
    return mtar(m,h,om,w)*(1+(1+w)/(1-3*w)*(1/OmegaM(1,om,w)-1))**(1/3)

#Função alfa(M) proposta
def alfa(m):
    return 0.05*(1 - 0.2*(np.log10(m/1e11)))

def delta_c_m(m):
    return 1.686*(1+6.54*alfa(m))

#Não consegui entender bem como extrair os valores de alfa...

#Função alfa que obtive
def alfa_m(m):
    return 0.829*(m)**(-0.111)

#mass_range = np.geomspace(1e11,1e15,3)
#print(mass_range)
#print(alfa(mass_range))
#print(delta_c_m(mass_range))

#Função mtar com shear e rotation
def mtar_sr(m,h,om,w):
    return mtar(m,h,om,w)*(1 - alfa(m)*delta_c_m(m))**(1/3)

def diff(m,w):
    return 1 - mtar_de(m,1,0.315,w)/mtar(m,0.673,0.315,w)


#Função de Massa: Usaremos dados tabelados.
#Criar uma lista chamável M = [ LISTA ]
#Primeira opção: 
#dados = np.genfromtxt("Masstab.txt",dtype=("U16", float, float,float))
#print(dados[0])
#exit()
m_s  = np.loadtxt("Masstab.txt", dtype=float, skiprows=1, usecols=(1))
mass = np.loadtxt("Masstab.txt", dtype=float, skiprows=1, usecols=(2))
merP = np.loadtxt("Masstab.txt", dtype=float, skiprows=1, usecols=(3))
merM = np.loadtxt("Masstab.txt", dtype=float, skiprows=1, usecols=(4))
radi = np.loadtxt("Masstab.txt", dtype=float, skiprows=1, usecols=(5))
rerP = np.loadtxt("Masstab.txt", dtype=float, skiprows=1, usecols=(6))
rerM = np.loadtxt("Masstab.txt", dtype=float, skiprows=1, usecols=(7))

for i in range(len(mass)):
    merP[i]= m_s[i]*merP[i]
    merM[i]= m_s[i]*merM[i]
    mass[i]= m_s[i]*mass[i]
    
    print("{:4.4e}".format(mass[i]),"{:4.4e}".format(merP[i]),"{:4.4e}".format(merM[i]))



#print("{:4.4e}".format(mass[-2]))
#w_range=np.linspace(-4,-0.5,100)
#plt.plot(w_range, (mtar(mass[-2],0.673,0.315,w_range)))
#plt.show()
#plt.errorbar(radi, mass, xerr=rerP, yerr=merP, fmt='.k', capsize=3) 

mass_range = np.geomspace(1e11,5e15,100)

#w_range = [-0.6,-0.8,-1.0,-1.2  ,-1.4]
#color   = ["g" ,"k" ,"b" ,"pink","r"]

w_range = [-0.5,-1.0,-2.5]
color   = ["g" ,"b" ,"r"]


#plt.show()

#Estou implementando uma rotina semalhante a esta para colocar os erros assimétricos:
# link: https://matplotlib.org/stable/gallery/statistics/errorbar_features.html#sphx-glr-gallery-statistics-errorbar-features-py

M_asymmetric_error = [merM, merP]#Minha ideia aqui era associar o contador ic a lista.É possível?
R_asymmetric_error = [rerM, rerP]
x = mass
y = radi
#plt.errorbar(x,y,  xerr=0.1*x, yerr=0.1, fmt = '.', markersize = 4, markeredgewidth = 2, markerfacecolor ='k', capsize=5)

plt.errorbar(x,y, xerr = M_asymmetric_error, yerr = R_asymmetric_error, fmt = '.', markersize = 2, markeredgewidth = 2, markerfacecolor ='k',            capsize=2)

ic =0
for wx in w_range:
    plt.plot(mass_range, mtar(mass_range,0.673,0.315,wx),c=color[ic])
    plt.plot(mass_range, mtar_de(mass_range,0.673,0.315,wx),c=color[ic],linestyle="--")
    #plt.plot(mass_range, mtar_sr(mass_range,0.673,0.315,wx),c=color[ic], linestyle="-.")
    #plt.plot(mass_range, mtar_popolo(mass_range,0.673,0.315,wx),c=color[ic], linestyle="-.")
    ic+=1

plt.xscale("log", nonposx='clip')
plt.yscale("log", nonposy='clip')
plt.ylabel("R")
plt.xlabel("M")
plt.ylim(0.08,10)
plt.xlim(1e11,1e16)
plt.show()

ik = 0
for wp in w_range:
     #plt.plot(mass_range, mtar_sr(mass_range,0.673,0.315,wp),c=color[ik], linestyle="-.")
     ik+=1

plt.xscale("log", nonposx='clip')
plt.yscale("log", nonposy='clip')
plt.ylabel("R")
plt.xlabel("M")
#plt.ylim(1e22,1e26)
plt.xlim(1e11,1e16)
#plt.show()
exit()





plt.xscale("log",nonposx='clip')
plt.yscale("log",nonposy='clip')
plt.ylabel("R[Mpc]")
plt.xlabel("M[Msun]")
plt.ylim(0.08,15)
plt.xlim(1e11,1e16)

plt.subplot(222)
plt.plot(w_range, diff(mass_range, w_range), 'k')
plt.xscale("log", nonposx='clip')
plt.yscale("log", nonposy='clip')
plt.ylabel("R_de/R_m")
plt.xlabel("w")

plt.show()
#print(M_asymmetric_error, R_asymmetric_erro
exit()


with open('Masstab.txt','r') as arquivo:
	Masstab = arquivo.readlines()
#print (Masstab)


#Segunda opção: Pandas
Mass_tab = pd.DataFrame({"Name": ["M31/MW", "M81","IC323","NGC253","CenA/M83","Local group", "Fornax-Eridanus", "Virgo"], 
                          "Mass": [(2.4e12),(8.9e11),(2.0e11),(2.34e11),(2.00e12),(3.14e12),(1.92e14),(1.67e15)], 
                          "Sup_error": [1.0,2.2,1.2,0.99,0.42,0.57,0.00,0.31],
                           "Inf_error": [-1.0,-2.2,-1.0,-1.34,-0.33,-0.54,0.00, -0.35]})
#print(Mass_tab)


exit()
#Função R(a): Falta Mi que ainda não sei determinar ou se vamos usar dados tabelados'
def eqs(y,a,om,w):
    R,R_p = y
    dyda  = [R_p, R_p*((1/2*a)*(OmegaM(a,om,w) + OmegaE(a,om,w)*(1+3*w))) -
             OmegaM(a,om,w)*(R/2*a**2)*(3*Mi/(4*math.pi*densM(a,om,w)*R**3) + (1/OmegaM(a,om,w) -1)*(1+3*w))]
    #return dyda
    
#Condições iniciais:
y0 =[1,1]#Colocar R = 0? e R_p = 0?
ai = 0.001 # Ou devo parametrizar tudo com a0 = 1 e supor que em tal época as quantidades era uma fração do que são hoje?
a_grid = np.linspace(ai,1,1000) #Aqui já estou usando 0.001<= a <= 1'


#Solução via odeint:

def sol(om,w):
    sol = odeint(eqs, y0, a_grid, args=(om,w), rtol=1e-16, atol=1e-10, printmessg=0)
    R     = interp1d(a_grid, sol[:,0], kind='cubic')
    R_p   = interp1d(a_grid, sol[:,1], kind='cubic')
    return R,R_p
R1, R_p1 = sol(0.3,-1)   #om = 0.3 e w = -1 (Lambda)
R2, R_p2 = sol(0.34,-0.8)#om =0.34 e w = -0.8 (dark  energy)


Rt1 = TAR(0.3,-1)
Rt2 = TAR(0.34,-0.8)
Rt3 = TAR(2.5,-1.5)

#Plots do TAR
