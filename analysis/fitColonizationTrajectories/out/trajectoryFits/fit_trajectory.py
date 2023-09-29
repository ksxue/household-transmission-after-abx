import pylab
import numpy
from scipy.optimize import minimize, brentq
import sys

# Theory line for logistic model
def logistic_logfs(ts,r,K,tstar):

	first_term = numpy.log(K)
	second_term = (-r*(tstar-ts)-numpy.log1p(numpy.exp(-r*(tstar-ts)*(ts<tstar))))*(ts<tstar)

	third_term = (-1*numpy.log1p(numpy.exp(-r*(ts-tstar)*(ts>=tstar))))*(ts>=tstar)
	return first_term+second_term+third_term

# Theory line for piecewise model
def piecwise_logfs(ts,r,K,tstar):
	#print("Theory: r=%g, K=%g, tstar=%g" % (r,K,tstar))
	trajectory = numpy.log(K)+r*(ts-tstar)*(ts<tstar)
	#print(numpy.exp(trajectory))
	return trajectory

# Fitting function for logistic trajectory
def fit_logistic_trajectory(ts, fs, fmin=1e-05):

	
	# Loss function is squared deviation from theory line, but putting a cap at fmin
	loss_function = lambda x: numpy.square(numpy.log(fs)-numpy.fmax(logistic_logfs(ts,numpy.exp(x[0]),numpy.exp(x[1]),x[2]),numpy.log(fmin))).sum()+x[0]+0.1*x[1]
	
	# Need to make an initial guess
	lr0 = -1
	lK0 = numpy.log(fs.max())
	tstar0 = ts[numpy.nonzero(fs>1e-04)[0][0]]
	#print("Guess:", lr0,lK0,tstar0)
	
	# Do the minimization
	res = minimize(loss_function,numpy.array([lr0,lK0,tstar0]))
	#print(res.x,res.success)
	r = numpy.exp(res.x[0])
	K = numpy.exp(res.x[1])
	tstar = res.x[2]
	#print(r,K,tstar)
	
	return r,K,tstar
	
# Same thing but for piecwise
def fit_piecwise_trajectory(ts, fs, fmins=[],tstar_guess=-1):

	if len(fmins)==0:
		fmins = numpy.ones_like(fs)*1e-05
		
	print("Fitting")
	print(ts)
	print(fs)

	# first truncate
	
	# plus version
	#loss_function = lambda x: numpy.square(numpy.log(fs)-numpy.log(numpy.exp(theory_logfs(ts,numpy.exp(x[0]),numpy.exp(x[1]),x[2]))+fmin)).sum()+x[0]
	
	# min version
	# Add a prior distribution
	# epsilon(log r - logrmin)^2
	loss_function = lambda x: numpy.square(numpy.log(fs)-numpy.fmax(piecwise_logfs(ts,numpy.exp(x[0]),numpy.exp(x[1]),x[2]),numpy.log(fmins))).sum()+0.1*x[0]+0.03*numpy.square(x[1]/11+1)
	
	
	
	lr0 = -1
	lK0 = numpy.log(fs.max())
	if tstar_guess<0:
		tstar0 = ts[numpy.nonzero(fs>=1e-03)[0][0]]
	else:
		tstar0 = tstar_guess
		
	print("Initial guess:", numpy.exp(lr0),numpy.exp(lK0),tstar0)
	#print("Guess:", lr0,lK0,tstar0)
	res = minimize(loss_function,numpy.array([lr0,lK0,tstar0]))
	#print(res.x)
	lr = res.x[0]
	lK = res.x[1]
	tstar = res.x[2]
	r = numpy.exp(lr)
	K = numpy.exp(lK)
	
	#print(r,K,tstar)
	
	# In many cases, we lack the temporal sampling to fully constrain r.
	# So now let's calculate the largest r value that is still consistent with the data

	# Maximum reasonable r (grows from detection threshold to 10% in one day)
	lrmax = numpy.log(numpy.log(1e-01/1e-05)/1)

	# Freeze K and tstar, and just look as function of r
	loss_function_without_prior = lambda x: numpy.square(numpy.log(fs)-numpy.fmax(piecwise_logfs(ts,numpy.exp(x),K,tstar),numpy.log(fmins))).sum()
	
	loss_fit = loss_function_without_prior(lr)
	loss_max = loss_function_without_prior(lrmax)
	
	dloss_max = numpy.log(2)
	
	if loss_max < loss_fit + dloss_max:
		# Max r is still within error tolerance threshold...
		rmax = numpy.exp(lrmax)
	else:
		# Find the root
		lrstar = brentq(lambda x: loss_function_without_prior(x)-(loss_fit+dloss_max), lr,lrmax)
		rmax = numpy.exp(lrstar)
	
	return r,K,tstar,rmax
	

# takes the second half of the trajectory near the minimum
def truncate_trajectory(ts,fs):

	fmin = fs.min()
	min_idxs = numpy.nonzero((fs<=(fmin+1e-09)))[0]
	idx = min_idxs[len(min_idxs)//2]
	
	return ts[idx:],fs[idx:]

# Calculate idxs of times where trajectory goes from above fupper to below flower
# (reports first idx below flower)	
def calculate_drops(ts,fs,flower=3e-05,fupper=1e-03,min_t=0):
	
	above_idxs = numpy.nonzero((fs>fupper)*(ts>=min_t))[0]
	
	if len(above_idxs)==0:
		# No timepoints above!
		return numpy.array([]),numpy.array([])
	
	drops = []
	jumps = []
	for dummy_idx in range(0,len(ts)):
		
		# Get first timepoint above
		first_above_idx = above_idxs[0]
		jumps.append(first_above_idx)
		
		# Get first timepoint below after the first timepoint above
		below_idxs = numpy.nonzero((fs<flower)*(ts>=ts[first_above_idx]))[0]
	
		if len(below_idxs)==0:
			# No timepoints below after the first timepoint above
			return numpy.array(drops), numpy.array(jumps)
		
		# There has been a drop below, 
		# -> increment the counter
		first_below_idx = below_idxs[0]
		
		drops.append(first_below_idx)
		
		# Calculate next timepoint above
		above_idxs = numpy.nonzero((fs>fupper)*(ts>=ts[first_below_idx]))[0]
		
		if len(above_idxs)==0:
			# Never get back above
			return numpy.array(drops), numpy.array(jumps)
			
		
if __name__=='__main__':
	
	example_ts = numpy.array([33,34,35,36,37,38,39,44,50,57,64,267,268])*1.0
	example_log10fs = numpy.array([-5,-5,-5,-5,-5,-5,-5, -3, -1.86, -1.84, -2.23, -2.2, -2.18])
	example_fs = numpy.power(10,example_log10fs)
	
	example_ts = numpy.array([33,34,35,36,37,38,39,44])
	example_fs = numpy.array([1e-05,1e-05,1e-05,1e-05,1e-05,1e-05,1e-05,1e-02])
	#truncated_ts, truncated_fs = truncate_trajectory(example_ts, example_fs)
	#print(truncated_ts)
	
	pylab.semilogy(example_ts,example_fs,'k.-')
	
	
	r,K,tstar = fit_logistic_trajectory(example_ts, example_fs)
	theory_ts = numpy.linspace(example_ts[0],70,100)
	theory_fs = numpy.exp(logistic_logfs(theory_ts,r,K,tstar))
	pylab.plot(theory_ts,theory_fs,'r-')
	
	r,K,tstar = fit_piecwise_trajectory(example_ts, example_fs)
	theory_ts = numpy.linspace(example_ts[0],70,100)
	theory_fs = numpy.exp(piecwise_logfs(theory_ts,r,K,tstar))
	pylab.plot(theory_ts,theory_fs,'b-')
	
	
	#print(theory_fs[0])
	
	pylab.xlim([30,70])
	pylab.ylim([1e-06,1])
	pylab.savefig('output.pdf',bbox_inches='tight')
	