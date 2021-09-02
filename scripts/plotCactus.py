import pylab as pl

colors=['green','blue','red','brown','orange']
mkrs= ['o','^','+','x','D']
for config in [[['DPSampler',bayestimes[1]],['WAPS',bayestimes[0]]],'bayes'], [[['DPSampler',pstimes[1]],['WAPS',pstimes[0]]],'ps']:
	colo = 0
	for alg in config[0]: #['DPSampler w/o JT',bayestimes[2]]]:
		totalNumCompleted = 0
		timelist = {}
		for time in alg[1]:
			if time >= 1000:
				timelist[1000] = 1
			elif time in timelist:
				timelist[time] += 1
				totalNumCompleted += 1
			else:
				timelist[time] = 1
				totalNumCompleted += 1
		alltimes = timelist.keys()
		alltimes.sort()
		numCompletedBenchsForEachTime = []
		rt = 0
		for time in alltimes:
			rt = rt + timelist[time]
			numCompletedBenchsForEachTime = numCompletedBenchsForEachTime + [rt]
		print('Total Completed Benchmarks for ',alg[0],' is ',totalNumCompleted,numCompletedBenchsForEachTime[-2])
		assert(totalNumCompleted == numCompletedBenchsForEachTime[-2])	
		pl.plot(numCompletedBenchsForEachTime,alltimes,label=alg[0], linewidth=3.0, color=colors[colo])		
		colo = colo + 1		
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42	
	pl.legend(loc='lower right')
	pl.tick_params(labelsize=16)
	pl.ylim(0,1010)
	pl.xlim(0,1090)
	pl.xlabel('Number of benchmarks solved', fontsize=20)
	pl.ylabel('Time (in seconds)', fontsize=20)
	#pl.title('Comparison of Number of benchmarks solved', fontsize=24)
	pl.axhline(y=7200, color='r', linestyle='--')
	pl.savefig('graphs/cactus_completed_'+config[1]+'.eps',bbox_inches='tight')
	pl.show()