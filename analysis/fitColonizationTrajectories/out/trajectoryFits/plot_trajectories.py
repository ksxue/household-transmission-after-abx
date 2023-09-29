import fit_trajectory
import parse_data
import pylab
import numpy
import sys

N = 1e13

# Denotes values of fit params to use when fitting fails
fit_params_failure = (0,0,-1)
fitted_parameters = {}

def calculate_last_preabx_timepoint(ts):
	last_preabx_idx = numpy.nonzero(ts < parse_data.FIRST_ABX_TIMEPOINT)[0][-1]
	return ts[last_preabx_idx]	

katherine_data = parse_data.parse_katherine_data()

desired_name="XIC-Bacteroides_nordii_55557"

# Data structure that holds trajectories for each subject/species combination
trajectories = {}
num_total_pops = 0
for subject in sorted(katherine_data):
	for species in sorted(katherine_data[subject]):
		num_total_pops += 1
		name = "%s-%s" % (subject,species)
		
		ts = katherine_data[subject][species][0]
		fs = katherine_data[subject][species][1]
		fmins = katherine_data[subject][species][2]
		
		last_abx_idx = numpy.nonzero(ts <= parse_data.LAST_ABX_TIMEPOINT)[0][-1]
		first_postabx_idx = numpy.nonzero(ts > parse_data.LAST_ABX_TIMEPOINT)[0][0]
		
		first_postabx_timepoint = ts[first_postabx_idx]
		
		last_preabx_idx = numpy.nonzero(ts < parse_data.FIRST_ABX_TIMEPOINT)[0][-1]
		first_abx_idx = numpy.nonzero(ts >= parse_data.FIRST_ABX_TIMEPOINT)[0][0]
		# Now apply a bunch of filters for trajectories we don't want to analyze
		# (this is a quick and dirty set, so can add some more)
		
		pre_abx_frequency = numpy.median(fs[:last_preabx_idx+1])
		median_nonzero_frequency = numpy.median(fs[fs>(4*fmins)])
		
		#name = "%s-%g" % (name,median_nonzero_frequency)
		
		
		# IF THERE ARE TIMEPOINTS BEFORE ABX THAT ARE MORE THAN 10-fold below median-nonzero frequency.
		
		# Locations where the trajectory drops from above 1e-03 to below 1e-05...
		drops_total, jumps_total = fit_trajectory.calculate_drops(ts,fs)
		drops_postabx, jumps_postabx = fit_trajectory.calculate_drops(ts,fs,min_t=ts[last_abx_idx])
		
		DISRUPTED = (pre_abx_frequency > 3e-04) and (((drops_total<=first_postabx_idx)*(drops_total>=last_preabx_idx)).sum() > 0)
	
		COLONIZED = (pre_abx_frequency < 1e-05) and ((jumps_total>last_preabx_idx).sum() > 0) 
		
		# OLD VERSION
		#COLONIZED = (pre_abx_frequency < 1e-05) and (((fs>=1e-03)*(ts>ts[last_preabx_idx])).sum() > 0)
		
		num_drops_total = len(drops_total)
		num_drops_postabx = len(drops_postabx)
		
		#print(name, num_drops_total, num_drops_postabx)
		if ((fs>=1e-03)*(ts>ts[last_abx_idx])).sum() < 1:
			# species not at high enough abundance after ABX so skip fitting!
			#if num_drops_postabx==0:
			#	print("Wouldn't have dropped!")
			#	print(ts)
			#	print(fs)
			print(name, "failed for never having high abundance after ABX")
			print(ts)
			print(fs)
			r,K,tstar = fit_params_failure
			fitted_parameters[name] = (r,K,tstar,r, "NO_HIGH_ABUNDANCE_POST_ABX")
			continue
			
		#abx_idxs = (ts>=parse_data.FIRST_ABX_TIMEPOINT)*(ts<=parse_data.LAST_ABX_TIMEPOINT)
		#if abx_idxs.sum() == 0:
			# Not measured during ABX, so skip
		#	print(name, "failed for not having any ABX timepoints")
		#	continue
		
		if not (DISRUPTED or COLONIZED):
			print(name, "failed for not being disrupted or colonized")
			print(ts)
			print(fs)
			r,K,tstar = fit_params_failure
			fitted_parameters[name] = (r,K,tstar,r, "NOT_DISRUPTED_OR_COLONIZED")
			continue
		
			
		#if (fs[last_abx_idx] > 1e-04) and (fs[first_postabx_idx] > 1e-04):
		#	print(name, "failed for not having low abundance at end of ABX", ts[last_abx_idx],fs[last_abx_idx]) 
			# didn't really die out, so skip
		#	if not (DISRUPTED or COLONIZED):
		#		pass
		#	else:
		#		print("WOULD NOT have failed for not being colonized or disrupted")
		#		print("DISRUPTED=", DISRUPTED, "COLONIZED=", COLONIZED)
		#		print(ts)
		#		print(fs)
				
		#	continue
		
		#else:
		#	if not (DISRUPTED or COLONIZED):
		#		print(name, "WOULD HAVE failed for not being colonized or disrupted")
				
		#		print("DISRUPTED=", DISRUPTED, "COLONIZED=", COLONIZED)
		#		print(pre_abx_frequency, drops_total)
		#		print(ts)
		#		print(fs)
				
		#if num_drops_postabx>0:
		#	print(ts)
		#	print(fs)
			
		#print("PASSED!")
		# Passed all the filters, so add it to the list
		
		#if name!=desired_name:
		#	continue
			
		trajectories[name] = (ts,fs,fmins,drops_total, jumps_total)

		
print("Processing %d species (of %d)..." % (len(trajectories),num_total_pops))		
#sys.exit(0)

tstars = []
rs = []
Ks = []


for name in sorted(trajectories):
	
	print(name)
		
	ts,fs,fmins,drop_idxs,jump_idxs = trajectories[name]
	
	max_t = ts[-1]+1
	
	last_preabx_idx = numpy.nonzero(ts < parse_data.FIRST_ABX_TIMEPOINT)[0][-1]	
	first_postabx_idx = numpy.nonzero(ts > parse_data.LAST_ABX_TIMEPOINT)[0][0]	
	
	if len(drop_idxs)>0:
		drop_ts = ts[drop_idxs]
		if (drop_ts>ts[first_postabx_idx]).sum()>0:
			#print("Switching max_t")
			max_t = drop_ts[(drop_ts>ts[first_postabx_idx])][0]	
		else:
			pass
			#print("Not switching max_t", drop_idxs)
			
	if fs[ts>=parse_data.LAST_ABX_TIMEPOINT][0]<1e-04:
		# dead at end of abx
		min_t = parse_data.LAST_ABX_TIMEPOINT
		#print("Not there at end of ABX")
	elif numpy.median(fs[:last_preabx_idx+1])<1e-04:
		# a colonizer not there at beginning of ABX
		min_t = ts[last_preabx_idx]
		#print("Not there ate beginning")
	else:
		if len(drop_idxs)==0:
			print("PROBLEM: NO DROP IDXS FOR DISRUPTED GUY")
			print(ts)
			print(fs)
		# disrupted guy, use the first drop after pre_abx_idx
		#print("Disrupted, use the first drop after pre_abx_idx")
		min_t = ts[drop_idxs[drop_idxs>=last_preabx_idx][0]]
		
	# Truncate trajectory so we just look at the post abx part
	# but before the species goes extinct again
	good_idxs = (ts>=min_t)*(ts<max_t)
	
	if good_idxs.sum()==0:
		print("Error: no good timepoints", name)
		print(min_t,max_t)
		print(ts)
		print(fs)
		print(drop_idxs)
		# Should not happen. Exit so we can debug
		# (you can assume that this didn't happen when I ran it)
		sys.exit(1)
	
	# Truncate trajectory to the region we care about	
	truncated_ts = ts[good_idxs]
	truncated_fs = fs[good_idxs]
	truncated_fmins = fmins[good_idxs]
	
	# Since we're going to be fitting log(f), add a minimum freq so log doesn't go crazy
	truncated_fs = numpy.fmax(truncated_fs,truncated_fmins)
	
	# Fit the piecewise model [log f = log K - r(tstar-t) * step_function(tstar-t)]
	r,K,tstar,rmax = fit_trajectory.fit_piecwise_trajectory(truncated_ts, truncated_fs,fmins=truncated_fmins)
	
	if r==0 or K==0 or (tstar<0) or (tstar>2e03):
		# Fit is weird
		# record standard N/A params
		r=0
		K=0
		tstar = -1
		reason = "FITTING_FAILED"
	else:
		reason = "SUCCESS"
		
	if tstar < 0:
		# Sometimes it helps to try one more time with a better guess for tstar
		
		# Get first timepoint above threshold
		tstar_guess = truncated_ts[(truncated_fs>truncated_fmins)][0]
		print("Refitting with tstar0=", tstar_guess)
		r,K,tstar,rmax = fit_trajectory.fit_piecwise_trajectory(truncated_ts, truncated_fs,fmins=truncated_fmins,tstar_guess=tstar_guess)
		if r==0 or K==0 or (tstar<0) or (tstar>2e03):
			# Fit is STILL weird
			# record standard N/A params
			print("Refit is still weird!")
			r=0
			K=0
			tstar = -1
			reason= "FITTING_FAILED"
		else:
			reason= "SUCCESS"
	
	# Add output to list
	rs.append(r)
	tstars.append(tstar)
	Ks.append(K)
	fitted_parameters[name] = (r,K,tstar,rmax,reason)
	print(name,r,K,tstar)
	
	# Make a little plot for this species
	
	# Calculate inferred establishment time  
	if tstar<ts[first_postabx_idx]:
		tau = truncated_ts[0]
	else:
		tau = tstar - numpy.log(N*r*K)/r
		tau = max([ts[first_postabx_idx],tau])
		
		
	# Timepoints in idx form
	t_idxs = numpy.arange(0,len(ts))
	
	# Use the real timepoints so we can make a theory trajectory in fake time
	theory_ts = ts
	
	# Make theory curve for piecwise model
	piecwise_fs = numpy.exp(fit_trajectory.piecwise_logfs(theory_ts,r,K,tstar))
	
	# Make a plot in scaled time
	if len(t_idxs)>50:
		width=10
	else:
		width=3.42
	pylab.figure(figsize=(width,2))
	
	pylab.semilogy(t_idxs,numpy.fmax(fs,1e-06),'k.')
	pylab.semilogy(t_idxs,piecwise_fs,'r-')
	
	pylab.xticks(t_idxs,[str(t) for t in ts],rotation='vertical',fontsize=6)
	pylab.xlabel('Timepoint')
	pylab.ylabel('Rel. abundance')
	pylab.ylim([1e-06,1])
	pylab.xlim([-1,len(t_idxs)])
	
	pylab.fill_between(t_idxs,numpy.ones_like(ts)*1e-06,numpy.fmax(fmins,1e-06),color='0.8')
	
	first_abx_idx = t_idxs[ts>=parse_data.FIRST_ABX_TIMEPOINT][0]-0.25
	last_abx_idx = t_idxs[ts<=parse_data.LAST_ABX_TIMEPOINT][-1]+0.25
	
	pylab.fill_between([first_abx_idx, last_abx_idx],[1e-06,1e-06],[1,1],color='r',alpha=0.5)
	
	max_tau_idx = t_idxs[ts>=tau][0]-0.5
	pylab.plot([max_tau_idx,max_tau_idx],[1e-06,1],'r:')
	
	pylab.savefig('output_scaled/%s.pdf' % name,bbox_inches='tight')
	pylab.close(pylab.gcf())
	
	# Make a plot in real time
	
	# Make theory curve for piecwise model
	theory_ts = numpy.linspace(theory_ts[0],theory_ts[-1],1000)
	piecwise_fs = numpy.exp(fit_trajectory.piecwise_logfs(theory_ts,r,K,tstar))
	max_piecwise_fs = numpy.exp(fit_trajectory.piecwise_logfs(theory_ts,rmax,K,tstar))
	pylab.figure(figsize=(3.42,2))
	pylab.semilogy(ts,numpy.fmax(fs,1e-06),'k.')
	pylab.semilogy(theory_ts,piecwise_fs,'r-')
	pylab.semilogy(theory_ts,max_piecwise_fs,'r-',alpha=0.5)
	
	pylab.semilogy([tau,tau],[1e-06,1],'r:',linewidth=0.5)
	pylab.xlabel('Timepoint (days)')
	pylab.ylabel('Rel. abundance')
	pylab.ylim([1e-06,1])
	pylab.fill_between(ts,numpy.ones_like(ts)*1e-06,numpy.fmax(fmins,1e-06),color='0.8')

	pylab.fill_between([parse_data.FIRST_ABX_TIMEPOINT, parse_data.LAST_ABX_TIMEPOINT],[1e-06,1e-06],[1,1],color='r',alpha=0.5)
	pylab.xlim([0,tstar+30])
	
	pylab.savefig('output_days/%s.pdf' % name,bbox_inches='tight')
	pylab.close(pylab.gcf())
	
	
	

# Now plot a summary across all species we looked at	
tstars = numpy.array(tstars)
rs = numpy.array(rs)
Ks = numpy.array(Ks)

good_idxs = (tstars>0)

tstars = tstars[good_idxs]
rs = rs[good_idxs]
Ks = Ks[good_idxs]
# total time it took to recover since end of ABX
dts = numpy.fmax(tstars-parse_data.LAST_ABX_TIMEPOINT,1)
taus = numpy.fmax(parse_data.LAST_ABX_TIMEPOINT,tstars - numpy.log(N*rs*Ks)/rs)

rtaus = rs*taus 

rdts = rtaus - rs*parse_data.LAST_ABX_TIMEPOINT

# r*ts-tau*abx

pylab.figure(figsize=(3.42,2))
pylab.loglog(dts,1/(rs*dts),'.')
pylab.xlabel('Total recovery time, $T_{sat}-T_{abx}$')
pylab.ylabel('$T_{growth}/(T_{sat}-T_{abx})$')
pylab.gca().set_ylim([1e-03,1])
pylab.savefig('output_summary/rt_vs_tstar.pdf',bbox_inches='tight')
	
pylab.figure(figsize=(3.42,2))
pylab.semilogx(tstars-taus,rdts,'.')
pylab.plot(tstars-taus,numpy.zeros_like(tstars),'k:',linewidth=0.25)
 

#pylab.xlabel('Saturation time, $T_{sat}$')
#pylab.ylabel('Establishment time, $\\tau$')
pylab.savefig('output_summary/tau_vs_tstar.pdf',bbox_inches='tight')	

pylab.figure(figsize=(3.42,2))
theory_ts = numpy.logspace(0,3,100)
survivals = numpy.array([((tstars-parse_data.LAST_ABX_TIMEPOINT)>=t).sum() for t in theory_ts])
pylab.semilogx(theory_ts,survivals,'k-')
pylab.semilogx(theory_ts,numpy.ones_like(theory_ts)*len(tstars),'k:')
pylab.xlabel('Tstar')
pylab.ylabel('Survival')
pylab.savefig('output_summary/tstar_survivals.pdf',bbox_inches='tight')	

output_file = open("output_summary/trajectory_fits.txt","w")
output_file.write("\t".join(["Subject","Species","r","K","tstar","Reason"]))

total_output = 0
good_output = 0
for name in sorted(fitted_parameters):
	r,K,tstar,rmax,reason = fitted_parameters[name]
	items = name.split("-")
	subject = items[0]
	species = items[1]
	
	# Don't output if a weird fit...
	if tstar<0:
		print(name, "failed fitting:", reason)
		#continue
	else:
		good_output+=1
	
	total_output+=1 
		
	output_file.write("\n")
	output_file.write("\t".join([subject,species,"%g" % r, "%g" % K,"%g" % tstar, "%s" % reason]))
	
output_file.close()

print("Done! %d trajectories fit (out of %d)" % (good_output, total_output))