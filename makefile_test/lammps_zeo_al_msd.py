#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import fort_analysis as fort
import helper as hlp


# In[2]:


fname = "/project/palmer/Ankur/lammps_zeolite/nonreax/nvt_woal_wmeth_May3/gmx_counterpart/dump_nve_twofs_unwrap.dcd"
groname = "/project/palmer/Ankur/lammps_zeolite/nonreax/nvt_woal_wmeth_May3/gmx_counterpart/test.gro"
ndxname="/project/palmer/Ankur/lammps_zeolite/nonreax/nvt_woal_wmeth_May3/gmx_counterpart/test.ndx"
funit = 500
#CV = []


# In[3]:


atomid = hlp.getnameandid_fromndx(ndxname)
#masses = hlp.getmass(groname)
#print masses
tatom, tframe, box = fort.readheaderdcd(fname)
print tatom, type(tatom)
print type(atomid['11']),atomid['11'].shape,len(atomid['11'])
print atomid['11']


# In[4]:


centers = atomid['11'][0:1]
#Hcenters1 = np.zeros((centers.shape),dtype=int)
#Hcenters1 = np.zeros((centers.shape),dtype=int)
#Hcenters1 = centers+1
#Hcenters2 = centers+2

#shell = atomid['OW']
#print shell
#Hshell1 = np.zeros((shell.shape),dtype=int)
#Hshell2 = np.zeros((shell.shape),dtype=int)
#Hshell1 = shell+1
#Hshell2 = shell+2


# In[ ]:


maxdtsteps = 500
dt_3dmsd = np.linspace(2,1000,500)
dr_3dmsd = np.zeros(dt_3dmsd.shape)
dx_3dmsd = np.zeros(dt_3dmsd.shape)
dy_3dmsd = np.zeros(dt_3dmsd.shape)
dz_3dmsd = np.zeros(dt_3dmsd.shape)
n_3dmsd = np.zeros(dt_3dmsd.shape)

for step in range(1,tframe-500):
    print ("Time: %.2f ps" % (step*0.5))
    dr_dt_3dmsd, dx_dt_3dmsd, dy_dt_3dmsd, dz_dt_3dmsd, n_dt_3dmsd = fort.msd(atomid['11'],fname,step,1,tframe,tatom,maxdtsteps,len(atomid['11']))
    dr_3dmsd = dr_3dmsd + dr_dt_3dmsd
    dx_3dmsd = dx_3dmsd + dx_dt_3dmsd
    dy_3dmsd = dy_3dmsd + dy_dt_3dmsd
    dz_3dmsd = dz_3dmsd + dz_dt_3dmsd

onebylen = 1.0/len(atomid['11'])


# In[ ]:


np.savetxt("lmp_zeo_woal.txt",(dr_3dmsd, dx_3dmsd, dy_3dmsd, dz_3dmsd, dt_3dmsd*0.5, n_3dmsd), header = "Angstrom^2 , ps , number")

