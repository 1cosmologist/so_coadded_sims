import numpy as np 

spectra_path = '/global/cfs/cdirs/sobs/delensing/official/sim_cls/'

Cl_unlensed = np.loadtxt(spectra_path+'cosmo2017_10K_acc3_scalCls.dat')
lmin_unlensed = int(np.min(Cl_unlensed[:,0]))
lmax_unlensed = int(np.max(Cl_unlensed[:,0]))

Cl_lensed = np.loadtxt(spectra_path+'cosmo2017_10K_acc3_lensedCls.dat')
lmin_lensed = int(np.min(Cl_lensed[:,0]))
lmax_lensed = int(np.max(Cl_lensed[:,0]))

r_ini = 1.
Cl_tensor = np.loadtxt('../resources/cosmo2017_10K_acc3_camb_tensorCls_r1.dat')
lmin_tensor = int(np.min(Cl_tensor[:,0]))
lmax_tensor = int(np.max(Cl_tensor[:,0]))

lmax_tensor = lmax_lensed
lmax_tot = np.max([lmax_lensed, lmax_tensor])

lensedDl = np.zeros((lmax_tot+1,4)) 
lensedDl[lmin_lensed:np.min([lmax_lensed, lmax_tot])+1,0:4] = Cl_lensed[:np.min([lmax_lensed, lmax_tot]) - lmin_lensed + 1,1:5]

unlensedDl = np.zeros((lmax_tot+1,4)) 
unlensedDl[lmin_lensed:np.min([lmax_unlensed, lmax_tot])+1,0:2] = Cl_unlensed[:np.min([lmax_unlensed, lmax_tot]) - lmin_unlensed + 1,1:3]
unlensedDl[lmin_lensed:np.min([lmax_unlensed, lmax_tot])+1, 3]  = Cl_unlensed[:np.min([lmax_unlensed, lmax_tot]) - lmin_unlensed + 1,3]

tensorDl = np.zeros((lmax_tot+1,4)) 
tensorDl[lmin_tensor:np.min([lmax_tensor, lmax_tot])+1,0:4] = Cl_tensor[:np.min([lmax_tensor, lmax_tot]) - lmin_tensor + 1,1:5]


def get_total(lmax=lmax_tot, r=1, Alens=1.):
    if lmax > lmax_tot:
        print(f"ERROR: lmax is too high! Set anything upto {lmax_tot}.")
        return None
    return (Alens * lensedDl + (r/r_ini) * tensorDl)[:lmax+1]

def get_lensed_scalar(lmax=lmax_tot, Alens=1.):
    if lmax > lmax_tot:
        print(f"ERROR: lmax is too high! Set anything upto {lmax_tot}.")
        return None
    return (Alens * lensedDl)[:lmax+1]

def get_tensor(lmax=lmax_tensor, r=1):
    if (lmax > lmax_tot):
        print(f"ERROR: lmax is too high! Set anything upto {lmax_tensor}.")
        return None
    if (lmax > lmax_tensor): print(f"WARNING: lmax is higher than lmax_tensor={lmax_tensor}. Padding with zeros!")

    return ((r/r_ini) * tensorDl)[:lmax+1]

def get_unlensed_scalar(lmax=lmax_tot):
    if lmax > lmax_tot:
        print(f"ERROR: lmax is too high! Set anything upto {lmax_tot}.")
        return None
    return unlensedDl[:lmax+1]



