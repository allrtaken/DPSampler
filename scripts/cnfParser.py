import sys

def parseCNF(inF):
	print('Parsing file ',inF)
	weighted = False
	projected = False
	clauses = []
	nVars = -1
	nCls = -1
	f = open(inF,'r')
	clsCounter = 0
	pEncountered = False
	wEncountered = False
	vpEncountered = False
	litWts = {}
	projVars = set()
	for line in f:
		if line.startswith('c'):
			continue
		if line.startswith('p'):
			if pEncountered:
				print ("Already encountered line starting with p. Second Line: ",line)
				sys.exit(1)
			else:
				pEncountered = True
			wds = line.split()
			if wds[1].startswith('w'):
				weighted = True
			if wds[2][1] == 'p':
				projected = True
			nVars = int(wds[2])
			nCls = int(wds[3])
			continue
		if line.startswith('w'):
			if weighted == False and wEncountered == False:
				print ("WARNING: weights encountered in a file with header p without wcnf. Ignoring weights..")
			wEncountered = True
			wds = line.split()
			litWts[int(wds[1])] = float(wds[2])
			assert(len(wds)<4 or wds[3] == '0')
			continue
		if line.startswith('vp'):
			if projected == False and vpEncountered == False:
				print ("WARNING: projection vars encountered in a file with header p without pcnf. Ignoring weights..")
				vpEncountered = True
				continue
			wds = line.split()
			assert(wds[-1] == '0')
			for i in range(1, len(wds)):
				projVars.add(int(wds[i]))
			continue
		if not pEncountered:
			print ("Encountered clause line:",line," before p cnf..")
			sys.exit(1)
		clsCounter += 1
		wds = line.split()
		cl = []
		for wd in wds:
			cl.append(int(wd))
		assert(cl[-1]==0)
		cl.pop()
		clauses.append(cl)
	f.close()
	print('File parsed!')
	return clauses,weighted,projected, litWts, projVars, nVars, nCls

