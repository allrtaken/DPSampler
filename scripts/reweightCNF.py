import sys, os, random

from cnfParser import parseCNF


inCNF = sys.argv[1]
numWeights = int(sys.argv[2])
outF = sys.argv[3]

probLow = 0.3
probHi = 0.7

clauses,weighted,projected, litWts, projVars, nVars, nCls = parseCNF(inCNF)

if weighted == True:
	print("WARNING: Input formula is already weighted. Ignoring those weights\n")

if projected:
	pVars = sorted(projVars)
else:
	pVars = range(1,nVars+1)

wtLits = random.sample([v*k for v in pVars for k in (-1,1)], numWeights)
wtLits.sort()

wts = []
for lit in wtLits:
	wts.append(random.uniform(probLow,probHi))

f = open(outF,'w')
if projected:
	f.write('p wpcnf '+str(nVars)+' '+str(nCls)+'\n')
	f.write('vp ')
	for v in pVars:
		f.write(str(v)+ ' ')
	f.write('0\n')
else:
	f.write('p wcnf '+str(nVars)+' '+str(nCls)+'\n')
for i in range(len(wts)):
	wt = wts[i]
	lit = wtLits[i]
	f.write('w '+str(lit)+' '+str(wt)+'\n')
for cl in clauses:
	for lit in cl:
		f.write(str(lit)+' ')
	f.write('0\n')
f.close()