import numpy as np

def readgro(getatomnames,gronames):
	grodat = []
	atomnames = getatomnames.split(',')
	print atomnames
	atomid = {}
	for i in atomnames:
		atomid[i] = [] 
	f = open(gronames.split()[0],'r')
	grodat = f.read().split('\n')[2:-2]
	f.close()
	for line in grodat:
		linesplit = line.split()
		if len(linesplit)==6 or len(linesplit)==9:
			if ('W' in linesplit[1]):
				linesplit[1] = linesplit[1][:2]
			try:
				atomid[linesplit[1]] += [int(linesplit[2])]
			except KeyError:
				continue
		else:
			id = linesplit[1][-5:]
			if ('W' in linesplit[1]):
				linesplit[1] = linesplit[1][:2]
			try:
				atomid[linesplit[1]] += [int(id)]
			except KeyError:
				continue
	for key in atomid:
		atomid[key] = np.array(atomid[key])
	return atomid


# In[ ]:


def getidfromanglesitp(itpfile):
	f = open(itpfile,'r')
	itpdata = f.read().split('\n')
	f.close()
	OH = [] #np.zeros(0,dtype=np.int32)
	HOH1 = [] #np.zeros(0,dtype=np.int32)
	HOH2 = [] #np.zeros(0,dtype=np.int32)
	for line in itpdata:
		linesplit = line.split()
		if len(linesplit) == 6 :
			OH += [int(linesplit[1])] 
			HOH1 += [int(linesplit[2])] 
			HOH2 += [-1] 
	return OH,HOH1,HOH2


def getidfromndxfile(ndxname):
	f = open(ndxname,'r')
	atomi = np.array(f.read().split()).astype(np.int32)
	f.close()
	atomid2 = []
	for i in atomi:
		atomid2 += [i]

def getnameandid_fromndx(ndxname):
	f = open(ndxname,'r')
	data = f.read().split('[')
	f.close()
	atomid = {}
	for dat in data:
		try:	
			dat2 = dat.split(']')
			atomid[dat2[0].strip().split('NAME_')[-1]] = np.array(dat2[1].strip().split()).astype(np.int32)
		except IndexError:
			continue
	return atomid

def writeCVhead(heads,fname):
	with open(fname,'w') as f:
		f.write('# '+' '.join(heads)+'\n')
def writeCVstyle(arr,fname):
	with open(fname,'a+') as f:
		f.write(' '+' '.join(str(val) for val in arr)+'\n')
def getmass(groname):
	mass = {}
	mass['OW'] = 16.
	mass['HW'] = 1.008
	mass['MW'] = 0.
	mass['Si'] = 28.08550
	mass['H'] = 1.00794
	mass['O'] = 15.9994
        mass['CH3'] = 15.034
        mass['CH2'] = 14.026
        mass['OH'] = 15.9994
        mass['HO'] = 1.008 
	f = open(groname.split()[0],'r')
	grodat = f.read().split('\n')
	totatom = int(grodat[1].strip())
	massarray = np.zeros((totatom))
	grodat = grodat[2:-2]
	f.close()
	for line in grodat:
		linesplit = line.split()
		if len(linesplit)==6 or len(linesplit)==9:
			if ('W' in linesplit[1]):
				linesplit[1] = linesplit[1][:2]
			massarray[int(linesplit[2])-1] = mass[linesplit[1]]
		else:
			id = linesplit[1][-5:]
			if ('W' in linesplit[1]):
				linesplit[1] = linesplit[1][:2]
			massarray[int(id)-1] = mass[linesplit[1]]
	return massarray
