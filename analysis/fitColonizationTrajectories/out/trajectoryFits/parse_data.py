import numpy

FIRST_ABX_TIMEPOINT = 29
LAST_ABX_TIMEPOINT = 34

def parse_katherine_data(filename="speciesRecoveryColonizationStrainTurnovers-trajectories.txt"):
	
	katherine_data = {}
	
	file = open(filename,"r")
	file.readline()
	for line in file:
		items = line.split()
		subject = items[1]
		t = int(items[3])
		species = items[4]
		f = float(items[5])
		
		# Estimate of minimum possible frequency (half of one read)
		fmin = 0.5/float(items[6]) 
		
		if subject not in katherine_data:
			katherine_data[subject] = {}
			
		if species not in katherine_data[subject]:
			katherine_data[subject][species] = {}
			
		katherine_data[subject][species][t] = (f,fmin)
	file.close()
		
	new_katherine_data = {}
	for subject in sorted(katherine_data):
		new_katherine_data[subject] = {}
		for species in sorted(katherine_data[subject]):
			
			ts = []
			fs = []
			fmins = []
			for t in sorted(katherine_data[subject][species]):
				f,fmin = katherine_data[subject][species][t]
				ts.append(t)
				fs.append(f)
				fmins.append(fmin)
				
			new_katherine_data[subject][species] = (numpy.array(ts),numpy.array(fs),numpy.array(fmins))
			
		
	return new_katherine_data