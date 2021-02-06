import numpy as np
import fort_analysis as fort
import helper as hlp

fname = "../npt.dcd"
groname = "../npt.gro"
ndxname="../npt.ndx"
funit = 500
CV = []
avgr = np.array([0.,0.,0.])
costheta = 0.

atomid = hlp.getnameandid_fromndx(ndxname)
masses = hlp.getmass(groname)
print masses
tatom, tframe, box = fort.readheaderdcd(fname)
msg = fort.opendcd(funit,fname)
msg = fort.skipheaderdcd(funit)
print tatom, type(tatom)
print type(atomid['OW']),atomid['OW'].shape
for step in range(tframe):
	rcoord,box2 = fort.readframedcd(funit,tatom)
	binv = 1.0/box2
	#print type(rcoord), rcoord.shape
	avgr = fort.average_loc(rcoord,
				atomid['OW'],
				16*np.ones(len(atomid['OW']),
				dtype=int),
				len(atomid['OW']),tatom)
	costheta = fort.endtoend_orient(box2,binv,
					rcoord[0,:],
					[1,2,3],rcoord,3,tatom)
	#print avgr,costheta
fort.closedcd(funit)
