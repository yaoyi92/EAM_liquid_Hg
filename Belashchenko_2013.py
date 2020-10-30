### The liquid mercury EAM potential follow Belashchenko 2013
### Belashchenko, David Kirillovich. "Application of the embedded atom model to liquid mercury." High Temperature 51.1 (2013): 40-48.
### Also note the http://towhee.sourceforge.net/forcefields/belash2013.html which point out some of the typos in the paper.
### units: Angstrom and eV

#from ase.calculators.eam import EAM
import numpy as np

#eam_Hg = EAM(potential='Al_zhou.eam.alloy')

## H(ri,ri+1) in eq.(3)
def H_fun(r,r1,r2):
    return np.heaviside(r-r1,0) - np.heaviside(r-r2,0)

## eq.(3)
def phi(r):
    output = 0.0
    n = 5
    k = 8
    r_partition = [2.55,3.10,3.70,4.60,6.20,8.30]
    a = [[-0.58730732649565e-1, -0.8725311607224e-1, 0.23383679109512e1, 0.54150014961337e2,
           0.50041303903954e3,  0.22115293102511e4, 0.50976873712314e4, 0.58975308220171e4, 0.27077049182679e4],
         [-0.55649321526289e-1, 0.18427524715662e-1, 0.41591385641041e0, 0.10555632657899e2, 
           0.94856510285549e2, 0.38745647613661e3, 0.80080059311747e3, 0.82197252535534e3, 0.33467726782330e3],
         [-0.14076511375606e-1, 0.44487200677395e-1, -0.85710953624001e-2, -0.17535348275128e0,
          -0.18290639200802e1, -0.61396850926303e1, -0.93832108328515e1, -0.68240837932190e1, -0.19167854338407e1],
         [-0.50217746756971e-2, -0.25080999359488e-1, -0.57924106076762e-1, -0.21114344856450e0, 
          -0.63905669655741e0, -0.11066271548233e1, -0.10219966911147e1, -0.46628922422428e0, -0.82702468028215e-1],
         [0.0e0, 0.0e0, -0.11675312540717e0, -0.42406122320448e0,
          -0.71663215147580e0, -0.64344074447359e0, -0.31451262110966e0, -0.78960959177628e-1, -0.79477007838712e-2]
           ]
    for i in [1,2,3,4,5]:
        for m in [0,1,2,3,4,5,6,7,8]:
            output += a[i-1][m]*(r-r_partition[i])**m*H_fun(r,r_partition[i-1],r_partition[i])
    output += (0.169356 - 6.34432*(2.55-r) + 3.8 * (np.exp(1.96*(2.55-r))-1))* (1 - np.heaviside(r-2.55,0))
    return output

## eq.(6)
def psi(r):
    p1 = 4.8019
    p2 = 1.3095
    return p1 * np.exp(-p2*r)


## eq.(7)-(11)
def phi_embed(rho):
    rho_partition = [1.0, 0.89, 0.81, 0.71, 0.62, 0.55, 0.47, 0.30, 1.20, 2.80]
    a = [-0.08798, -0.078461, -0.073575, -0.058668, -0.051551, -0.038358, -0.005728, -0.056069, -0.056512, 2.625839]
    b = [np.nan, -0.173074, 0.050926, -0.349074, 0.190926, -0.567874, -0.247874, 0.840126, 0.314680, 2.629723]
    c = [0.7867,-1.40,2.00,-3.00,5.42,-2.00,-3.20,4.00,0.980,0.230]
    m = 1.70
    n = 3.00
    output = 0.0
    output += (a[0] + c[0] * (rho - rho_partition[0])**2) * H_fun(rho,rho_partition[1], rho_partition[8])
    for i in [2,3,4,5,6,7]:
        output += (a[i-1] + b[i-1] * (rho - rho_partition[i-1]) + c[i-1] * (rho - rho_partition[i-1])**2) * \
                   H_fun(rho,rho_partition[i],rho_partition[i-1])
    output += (a[7] + b[7] * ( rho - rho_partition[7]) + c[7] * (rho - rho_partition[7])**2) * \
                   (2.0 * rho / rho_partition[7] - (rho / rho_partition[7])**2) * \
                   (1 - np.heaviside(rho - rho_partition[7],0))
    if (rho <= rho_partition[9]) and (rho > rho_partition[8]):
        output += a[8] + b[8] * (rho - rho_partition[8]) + c[8] * (rho - rho_partition[8])**m
    output += (a[9] + b[9] * (rho - rho_partition[9]) + c[9] * (rho - rho_partition[9])**n ) * \
                   np.heaviside(rho - rho_partition[9],0)
    return output


##### write the eam potential in eam.alloy (setfl) format
##### description of the format https://lammps.sandia.gov/doc/pair_eam.html
#####
print("###CONTRIBUTOER: Yi Yao yaoyi92 AATT gmail.com CITATION: Belashchenko, David Kirillovich.  High Temperature 51.1 (2013): 40-48.")
print("#")
print("#")
Nelement = 1
element = "Hg"
print(Nelement, element)
Nrho = 10001
drho = 0.002
Nr = 10001
dr = 0.001
cutoff = 8.5
print(Nrho, drho, Nr, dr, cutoff)
atomic_number = 80
mass = 200.59
lattice = 1.0
lattice_type = "FCC"
print(atomic_number, mass, lattice, lattice_type)
rho_list = np.linspace(0,(Nrho-1)*drho,Nrho)
for i_rho in range(Nrho):
    print(phi_embed(rho_list[i_rho]))
r_list = np.linspace(0,(Nr-1)*dr,Nr)
for i_r in range(Nr):
    print(psi(r_list[i_r]))
for i_r in range(Nr):
    print(r_list[i_r]*phi(r_list[i_r]))
#####


##### some test for debug

#rho_test = np.linspace(0.0,5.0,10001)
#phi_embed_test = []
#for rho in rho_test:
#    phi_embed_test.append(phi_embed(rho))
#for i in range(10001):
#    print(rho_test[i], phi_embed(rho_test[i]))

#r_test = np.linspace(0,9,9001)
#for i in range(9001):
#    print(r_test[i], phi(r_test[i]))
#phi_test = []
#for r in r_test:
#    phi_test.append(phi(r))
