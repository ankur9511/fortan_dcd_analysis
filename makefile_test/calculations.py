import numpy as np
import fort_analysis as fort
import helper as hlp

fname = "npt.dcd"
groname = "npt.gro"
ndxname="npt.ndx"
funit = 500
CV = []

atomid = hlp.getnameandid_fromndx(ndxname)
masses = hlp.getmass(groname)
#print masses
tatom, tframe, box = fort.readheaderdcd(fname)
msg = fort.opendcd(funit,fname)
msg = fort.skipheaderdcd(funit)
print tatom, type(tatom)
print type(atomid['OW']),atomid['OW'].shape,len(atomid['OW'])

centers = atomid['OW'][0:1]
Hcenters1 = np.zeros((centers.shape),dtype=int)
Hcenters1 = np.zeros((centers.shape),dtype=int)
Hcenters1 = centers+1
Hcenters2 = centers+2

shell = atomid['OW']
print shell
Hshell1 = np.zeros((shell.shape),dtype=int)
Hshell2 = np.zeros((shell.shape),dtype=int)
Hshell1 = shell+1
Hshell2 = shell+2

for step in range(tframe):
	rcoord,box2 = fort.readframedcd(funit,tatom)
	binv = 1.0/box2
	shellid,nhbpair,nshell = fort.nhbond_nshell(centers,Hcenters1,Hcenters2,
						  shell,Hshell1,Hshell2,
						  rcoord,box2,binv,
						  12.25,np.cos(30.0*np.pi/180.0),
						  tatom,len(centers),len(shell))
	print step,nhbpair
fort.closedcd(funit)
