import sys, os
from pysat.formula import CNF
from cnfParser import parseCNF
from getCNFSols import generateSols
from numpy.random import default_rng
import numpy as np
import math

inCNF = sys.argv[1]
outF = sys.argv[2]
nSamples = int(sys.argv[3])

logC = False

assert('TMP' in os.environ)
fname = os.path.basename(inCNF)

outF = os.environ['TMP']+'/'+fname+'.idealsamples'
print('Setting output file to ',outF)

def logsumexp(x):
	c = x.max()
	return np.sum(np.exp(x - c)), c



clauses,weighted,projected, litWts, projVars, nVars, nCls = parseCNF(inCNF)
formula = CNF()
formula.extend(clauses)

if projected:
	pVars = sorted(projVars)
else:
	pVars = range(1,formula.nv+1)

if not weighted:
	lWts = lambda x: 1
else:
	lWts = lambda x: litWts[x] if x in litWts else 1-litWts[-x] if -x in litWts else 1

num, sols, wts = generateSols(formula,pVars,lWts, logC)
if logC:
	totW,c = logsumexp(np.array(wts))
	print("\nTotal solutions:",num,' total weight:',totW*np.exp(c))
	probs = [np.exp(wt-c)/totW for wt in wts]
else:
	totW = sum(wts)
	print("\nTotal solutions:",num,' total weight:',totW)
	probs = [wt/float(totW) for wt in wts]

	
rng = default_rng()
samples = rng.choice(a=num,size=nSamples,p=probs)

of = open(outF,'w')
for i in range(nSamples):
	of.write(str(i+1)+', ')
	for j in sols[samples[i]]:
		of.write(str(j)+' ')
	of.write('\n')
of.close()