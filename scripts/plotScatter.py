import pylab as pl

colors=['green','blue','red','brown','orange']
mkrs= ['o','^','+','x','D']
for config in [bayestimes[0],bayestimes[1],bayesTWs,'Bayes'],[pstimes[0],pstimes[1],psTWs,'Pseudoweighted']:
	pl.scatter(config[2],config[0],label=alg[0], color='r', marker='o',facecolors='none')
	pl.scatter(config[2],config[1],label=alg[0], color='b', marker='x')
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42			
	pl.tick_params(labelsize=16)
	pl.ylabel('Time (sec)', fontsize=20)
	pl.xlabel('Tree-width', fontsize=20)
	pl.title('Time vs Treewidth for '+config[3], fontsize=24)
	pl.show()
for config in [bayestimes[0],bayestimes[1],'Bayes'],[pstimes[0],pstimes[1],'Pseudoweighted']:
	pl.scatter(config[1],config[0],label=alg[0], color='b', marker='x')
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42			
	pl.tick_params(labelsize=16)
	pl.xlabel('DPSampler Time (sec)', fontsize=20)
	pl.ylabel('WAPS Time (sec)', fontsize=20)
	pl.title('DPSampler Time vs Treewidth for '+config[2], fontsize=24)
	pl.show()