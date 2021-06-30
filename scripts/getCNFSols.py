from pysat.solvers import Solver
import pickle
import numpy as np

def generateSols(formula, pVars, litWts, logC):
	print('Generating solutions.. ',end='')
	sols = []
	wts = []
	num = 0
	with Solver(name='m22', bootstrap_with=formula) as oracle:
		for sol in oracle.enum_models():
			projSol = []
			if logC:
				wt = 0
			else:
				wt = 1
			for s in sol:
				if abs(s) in pVars:
					projSol.append(s)
					if logC:
						wt += np.log(litWts(s))
					else:
						wt *= litWts(s)
			sols.append(projSol)
			wts.append(wt)
			oracle.add_clause([-l for l in projSol])
			num += 1
			if num % 100 == 1:
				print (num, " ", end = '')
	return num,sols,wts

def writeSols(sols, outF):
	f = open(outF,'wb')
	pickle.dump(sols)
	f.close()

def readSols(inF):
	f = open(inF,'rb')
	sols = pickle.load(f)
	f.close()
	return sols