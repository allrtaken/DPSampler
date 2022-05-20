from sqlite3.dbapi2 import Cursor

from numpy.core.fromnumeric import mean
from processOutput import initDB, closeAll
import math, numpy as np
import pylab as pl
import sys

def geomean(xs):
	return math.exp(math.fsum(math.log(x) if x > 0 else math.log(x+0.3) for x in xs) / len(xs))

def geosd(xs, mu): #geometric standard deviation
	return math.exp(math.sqrt(math.fsum((math.log(x/mu))**2 if x>0 else (math.log((x+0.3)/mu))**2 for x in xs)/len(xs)))

def getStrings(tottime, title):
	wapsTimeStr = None
	dmcTimeStr = None
	titleStr = None
	if tottime == True:
		wapsTimeStr = 'waps.total_t'
		dmcTimeStr = 'dmc.total_t'
	else:
		wapsTimeStr = 'waps.compile_t'
		dmcTimeStr = '(dmc.totpreaddcomp_t + dmc.addcomp_t)'
	if title == 'all': #not including enc1.. do it separately if needed
		titleStr = '((waps.type=\'bayes\' and dmc.type=\'bayes\') OR (waps.type=\'ps\' and dmc.type=\'ps\') OR (waps.type=\'wb\' and dmc.type=\'wb\'))'
	elif title == 'allunique':
		titleStr = '((waps.type=\'bayes\' and dmc.type=\'bayes\') OR (waps.type=\'ps\' and dmc.type=\'ps\') OR (waps.type=\'wb\' and dmc.type=\'wb\'AND (waps.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')) ))'
	elif title=='bayes':
		titleStr = '(waps.type=\'bayes\' and dmc.type=\'bayes\')'
	elif title == 'ps':
		titleStr = '(waps.type=\'ps\' and dmc.type=\'ps\')'
	elif title == 'db': #dpmc benchmarks
		titleStr = '((waps.type=\'bayes\' and dmc.type=\'bayes\') OR (waps.type=\'ps\' and dmc.type=\'ps\'))'
	elif title == 'wb':
		titleStr = '(waps.type=\'wb\' and dmc.type=\'wb\')'
	else:
		print("Title not recognized. Exiting..")
		sys.exit(1)
	return wapsTimeStr,dmcTimeStr,titleStr
# (waps.type='bayes' and cudd.type='bayes') AND 

def calcSpeedups(cur:Cursor, tottime:bool, dmc_lib:str, title:str):
	wapsTimeStr,dmcTimeStr,titleStr = getStrings(tottime, title)
	rows = cur.execute("SELECT "+wapsTimeStr+"/"+dmcTimeStr+" from waps, "+dmc_lib+"_632ed69plus as dmc  WHERE "+titleStr+" AND (waps.name = dmc.name) AND (waps.type = dmc.type) AND ("+dmcTimeStr+" < 1000) AND ("+wapsTimeStr+" < 1000)").fetchall()
	mu = geomean(next(zip(*rows)))
	gsd = geosd(next(zip(*rows)), mu)
	print('geo-mean: ',mu,'  geo-sd:',gsd)

def calcFinished(cur:Cursor, tottime:bool, dmc_lib:str, title:str):
	wapsTimeStr,dmcTimeStr,titleStr = getStrings(tottime, title)
	both = cur.execute("SELECT COUNT(1) from waps, "+dmc_lib+"_632ed69plus as dmc  WHERE "+titleStr+" AND (waps.name = dmc.name) AND (waps.type = dmc.type) AND ("+dmcTimeStr+" < 1000) AND ("+wapsTimeStr+" < 1000)").fetchall()
	only_waps = cur.execute("SELECT COUNT(1) from waps, "+dmc_lib+"_632ed69plus as dmc  WHERE "+titleStr+" AND (waps.name = dmc.name) AND (waps.type = dmc.type) AND ("+dmcTimeStr+" IS NULL) AND ("+wapsTimeStr+" < 1000)").fetchall()
	only_dmc = cur.execute("SELECT COUNT(1) from waps, "+dmc_lib+"_632ed69plus as dmc  WHERE "+titleStr+" AND (waps.name = dmc.name) AND (waps.type = dmc.type) AND ("+dmcTimeStr+" < 1000) AND ("+wapsTimeStr+" IS NULL)").fetchall()
	print ("Both:",both[0][0]," only_waps:",only_waps[0][0]," only_dmc:",only_dmc[0][0])
	
def calcFastest(cur:Cursor, dmc_lib:str, title:str):
	# not computing for compile time because 'fastest' becomes blury. 
	# for eg we only consider ddnnf compilation time as compilation for waps, and not include time for annotation
	# as that is done in python (slow) and we want to be fair. so while its okay for num_finished, not much sense
	# for num_fastest 
	wapsTimeStr,dmcTimeStr,titleStr = getStrings(tottime=True, title=title)
	only_waps = cur.execute("SELECT COUNT(1) from waps, "+dmc_lib+"_632ed69plus as dmc  WHERE "+titleStr+" AND (waps.name = dmc.name) AND (waps.type = dmc.type) AND (( ("+dmcTimeStr+" < 1000) AND ("+wapsTimeStr+" < 1000) AND ("+dmcTimeStr+">"+wapsTimeStr+")) OR (("+wapsTimeStr+" < 1000) AND ("+dmcTimeStr+" IS NULL)))").fetchall()
	only_dmc = cur.execute("SELECT COUNT(1) from waps, "+dmc_lib+"_632ed69plus as dmc  WHERE "+titleStr+" AND (waps.name = dmc.name) AND (waps.type = dmc.type) AND (( ("+dmcTimeStr+" < 1000) AND ("+wapsTimeStr+" < 1000) AND ("+wapsTimeStr+">"+dmcTimeStr+")) OR (("+dmcTimeStr+" < 1000) AND ("+wapsTimeStr+" IS NULL)))").fetchall()
	print ("fastest_waps:",only_waps[0][0]," fastest_dmc:",only_dmc[0][0])

def plotSpeedupTWs(cur:Cursor, tottime:bool, dmc_lib:str, title:str):
	wapsTimeStr,dmcTimeStr,titleStr = getStrings(tottime, title)

	rows = cur.execute("SELECT "+wapsTimeStr+"/"+dmcTimeStr+", dmc.jt_width from waps, "+dmc_lib+"_632ed69plus as dmc  WHERE "+titleStr+" AND (waps.name = dmc.name) AND (waps.type = dmc.type) AND ("+dmcTimeStr+" < 1000) AND ("+wapsTimeStr+" < 1000) ORDER BY dmc.jt_width").fetchall()
	gms = []
	currJt = None
	currSUs = []
	for r in rows:
		newJt = r[1]
		if currJt == None: #first item
			currJt = newJt
			currSUs.append(r[0])
			continue
		if newJt != currJt:
			newGm = geomean(currSUs)
			gms.append((currJt,newGm, max(currSUs),min(currSUs)))
			currSUs = [r[0]]
			currJt = newJt
			continue
		assert newJt == currJt
		currSUs.append(r[0])
	colors=['green','blue','red','brown','orange']
	mkrs= ['o','^','+','x','D']
	pl.plot([g[0] for g in gms],[g[1] for g in gms],label='Geometric Means', linewidth=1.0, color=colors[0], marker='x')
	#pl.errorbar([g[0] for g in gms],[g[1] for g in gms],[[g[1]-g[3] for g in gms],[g[2]-g[1] for g in gms]])
	# for g in gms:
		# pl.errorbar(g[0],g[2]-g[3],bottom=g[3])
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42	
	pl.legend(loc='lower right')
	pl.tick_params(labelsize=16)
	#pl.ylim(0,1010)
	#pl.xlim(0,1090)
	pl.xlabel('Treewidths', fontsize=20)
	pl.ylabel('Speedup', fontsize=20)
	pl.title(title+' '+dmc_lib+' '+wapsTimeStr.split('.')[1], fontsize=24)
	pl.axhline(y=1, color='r', linestyle='--')
	#pl.savefig('graphs/cactus_completed_'+config[1]+'.eps',bbox_inches='tight')
	pl.show()

def processPar2Rows(rows):
	res = {}
	for r in rows:
		jtw = r[4]
		wTime = r[1] if r[1] != None else 2000
		cTime = r[2] if r[2] != None else 2000
		sTime = r[3] if r[3] != None else 2000
		name = r[0]
		if jtw in res:
			res[jtw][0].append(name)
			res[jtw][1].append(wTime)
			res[jtw][2].append(cTime)
			res[jtw][3].append(sTime)
		else:
			res[jtw] = [[name],[wTime],[cTime],[sTime]]
	return res

def plotPar2TWs(cur:Cursor, title:str, windowSize, avg): #, tottime:bool, dmc_lib:str, title:str):
	if title == 'all': #not including enc1.. do it separately if needed
		titleStr = '((waps.type=\'bayes\' and cudd.type=\'bayes\') OR (waps.type=\'ps\' and cudd.type=\'ps\') OR (waps.type=\'wb\' and cudd.type=\'wb\'))'
	elif title == 'allunique':
		titleStr = '((waps.type=\'bayes\' and cudd.type=\'bayes\') OR (waps.type=\'ps\' and cudd.type=\'ps\') OR (waps.type=\'wb\' and cudd.type=\'wb\'AND (waps.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')) ))'
	elif title=='bayes':
		titleStr = '(waps.type=\'bayes\' and cudd.type=\'bayes\')'
	elif title == 'ps':
		titleStr = '(waps.type=\'ps\' and cudd.type=\'ps\')'
	elif title == 'db':
		titleStr = '((waps.type=\'bayes\' and cudd.type=\'bayes\') OR (waps.type=\'ps\' and cudd.type=\'ps\'))'
	elif title == 'wb':
		titleStr = '(waps.type=\'wb\' and cudd.type=\'wb\')'
	else:
		print("Title not recognized. Exiting..")
		sys.exit(1)
	rows = cur.execute("SELECT waps.name, waps.total_t, cudd.total_t, sylvan.total_t, cudd.jt_width from waps, cudd_632ed69plus as cudd, sylvan_632ed69plus as sylvan  WHERE  (waps.name = cudd.name) AND (waps.name = sylvan.name) AND (waps.type = cudd.type) AND (waps.type = sylvan.type) AND "+titleStr+" AND cudd.jt_width > 0 ORDER BY cudd.jt_width").fetchall()
	timemap = processPar2Rows(rows)
	
	jtws = sorted(timemap.keys())
	# print(jtws,timemap[jtws[10]])
	l = len(jtws)
	assert l > windowSize and windowSize>0
	res = []
	for i in range(windowSize,l):
		res.append([])
		windows = [[],[],[],[]] #jtws, wtimes, ctimes, stimes
		for j in range(i-windowSize,i):
			windows[0].append(jtws[j])
			times = timemap[jtws[j]]
			windows[1] += times[1]
			windows[2] += times[2]
			windows[3] += times[3]
		res[-1] += [avg(windows[0]),avg(windows[1]),avg(windows[2]),avg(windows[3])]
		for i,val in enumerate(res[-1]):
			if math.isnan(val):
				print(i,windows,res[-1][i])
				sys.exit(1)
	# print(res)
	colors=['green','blue','red','brown','orange']
	mkrs= ['o','^','+','x','D']
	for item in [[1,'WAPS'],[2,'DPSampler']]:#,[3,'DPSampler w/ Sylvan']]:
	 	pl.plot([g[0] for g in res],[g[item[0]] for g in res],label=item[1], linewidth=1.0, color=colors[item[0]])

	# pl.plot([g[0] for g in res],[g[1] for g in res],label='WAPS', linewidth=1.0, color=colors[0])
	
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42	
	pl.legend(loc='lower right', fontsize=12)
	pl.tick_params(labelsize=14)
	#pl.ylim(0,1010)
	#pl.xlim(0,1090)
	pl.xlabel('Mean of '+str(windowSize)+' project-join treewidths',fontsize=16)#, fontsize=20)
	pl.ylabel('Mean PAR2 score of '+str(windowSize)+' widths',fontsize=16)#, fontsize=20)
	#pl.title(title+' '+dmc_lib+' '+wapsTimeStr.split('.')[1], fontsize=24)
	#pl.axhline(y=1, color='r', linestyle='--')
	pl.savefig('result_analysis/graphs/dmc_632ed69+/par2_'+title+'.eps',bbox_inches='tight')
	pl.show()

def plotCactus(cur:Cursor, title: str, caption:str, numBench): #, tottime:bool, dmc_lib:str, title:str):
	# assert title in ['bayes', 'ps']
	# not doing allunique because the graphs are plotted separately anyway
	if title=='all':
		title1 = 'bayes'
		title2 = 'ps'
		title3 = 'wb'
	elif title=='db':
		title1 = 'bayes'
		title2 = 'ps'
		title3 = 'ps'
	else:
		title1 = title
		title2 = title
		title3 = title
	rows_waps=cur.execute('SELECT t1.total_t, (SELECT COUNT(1) FROM waps t2 WHERE ((t2.type=\''+title1+'\') OR (t2.type=\''+title2+'\') OR (t2.type=\''+title3+'\'))  AND (t2.total_t <= t1.total_t)) AS num_complete FROM waps t1 WHERE ((t1.type=\''+title1+'\') OR (t1.type=\''+title2+'\') OR (t1.type=\''+title3+'\')) GROUP BY t1.total_t ORDER BY t1.total_t').fetchall()
	rows_cudd=cur.execute('SELECT t1.total_t, (SELECT COUNT(1) FROM cudd_632ed69plus t2 WHERE ((t2.type=\''+title1+'\') OR (t2.type=\''+title2+'\') OR (t2.type=\''+title3+'\')) AND (t2.total_t <= t1.total_t)) AS num_complete FROM cudd_632ed69plus t1 WHERE ((t1.type=\''+title1+'\') OR (t1.type=\''+title2+'\') OR (t1.type=\''+title3+'\')) GROUP BY t1.total_t ORDER BY t1.total_t').fetchall()
	# rows_sylvan=cur.execute('SELECT t1.total_t, (SELECT COUNT(1) FROM sylvan_632ed69plus t2 WHERE (t2.type=\''+title+'\') AND (t2.total_t <= t1.total_t)) AS num_complete FROM sylvan_632ed69plus t1 WHERE (t1.type=\''+title+'\') GROUP BY t1.total_t ORDER BY t1.total_t').fetchall()
	colors=['green','blue','red','brown','orange']
	mkrs= ['o','^','+','x','D']
	pl.plot([r[1] for r in rows_waps]+[rows_waps[-1][1]],[r[0] for r in rows_waps]+[1000],label='WAPS', linewidth=3.0, color=colors[2])
	pl.plot([r[1] for r in rows_cudd]+[rows_cudd[-1][1]],[r[0] for r in rows_cudd]+[1000],label='DPSampler', linewidth=3.0, color=colors[0])
	# pl.plot([r[1] for r in rows_sylvan],[r[0] for r in rows_sylvan],label='DPSampler w/ Sylvan', linewidth=3.0, color=colors[1])
	print([rows_waps[-1][1]], [rows_cudd[-1][1]])
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42	
	pl.legend(loc='upper left',fontsize=12)
	pl.tick_params(labelsize=16)
	pl.ylim(0,1010)
	# pl.yscale('log')
	pl.xlim(0,numBench+10)
	pl.ylabel('Total Time (seconds)', fontsize=20)
	pl.xlabel('Number of benchmarks solved', fontsize=20)
	#pl.title(caption+' Benchmark Set', fontsize=24)
	pl.axhline(y=1000, color='r', linestyle='--')
	pl.savefig('result_analysis/graphs/dmc_632ed69+/cactus_'+title+'.eps',bbox_inches='tight')
	pl.show()

def plotSamplingComparison():
	# see /home/adi/Downloads/dpsampler/results/dmc_632ed69+/sampling_time_comparison/notes.txt for data source
	data_x = [319,679,1004,1612,2075]
	data_y_num = [3.779, 30.97, 31.57, 265.36, 771.02]
	data_y_den = [0.211, 1.2, 2.046, 3.827, 5.452]
	data_y = [data_y_num[i]/data_y_den[i] for i in range(len(data_y_den))]
	colors=['green','blue','red','brown','orange']
	mkrs= ['o','^','+','x','D']
	pl.scatter(data_x,data_y, color=colors[0], marker='x')
	pl.plot(np.unique(data_x), np.poly1d(np.polyfit(data_x, data_y, 1))(np.unique(data_x)), label='Best-Fit Line')
	#pl.errorbar([g[0] for g in gms],[g[1] for g in gms],[[g[1]-g[3] for g in gms],[g[2]-g[1] for g in gms]])
	# for g in gms:
		# pl.errorbar(g[0],g[2]-g[3],bottom=g[3])
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42	
	pl.legend(loc='upper left',fontsize=12)
	pl.tick_params(labelsize=16)
	#pl.ylim(0,1010)
	#pl.xlim(0,1090)
	pl.xlabel('# nodes in Project-Join Tree', fontsize=20)
	pl.ylabel('Speedup', fontsize=20)
	#pl.title('Speedup of top-down sampling relative to bottom-up', fontsize=24)
	pl.axhline(y=1, color='r', linestyle='--')
	pl.savefig('result_analysis/graphs/dmc_632ed69+/sampling_comparison.eps',bbox_inches='tight')
	pl.show()

def plotSamplingComparison2(cur:Cursor, title: str):
	# bayes and ps comparison for iJCAI 22 (so no need for all unique)
	dmc1TimeStr = '(dmc1.sampgen_t)'
	dmc2TimeStr = '(dmc2.sampgen_t)'
	if title == 'all': #not including enc1.. do it separately if needed
		titleStr = '((dmc1.type=\'bayes\' and dmc2.type=\'bayes\') OR (dmc1.type=\'ps\' and dmc2.type=\'ps\'))'
	elif title=='bayes':
		titleStr = '(dmc1.type=\'bayes\' and dmc2.type=\'bayes\')'
	elif title == 'ps':
		titleStr = '(dmc1.type=\'ps\' and dmc2.type=\'ps\')'
	speedup_rows = cur.execute("SELECT "+dmc2TimeStr+"/"+dmc1TimeStr+" from cudd_632ed69plus dmc1, budmc dmc2  WHERE "+titleStr+" AND (dmc1.name = dmc2.name) AND (dmc1.type = dmc2.type) AND ("+dmc1TimeStr+" < 1000) AND ("+dmc2TimeStr+" < 1000)").fetchall()
	mu = geomean(next(zip(*speedup_rows)))
	gsd = geosd(next(zip(*speedup_rows)),mu)
	print ('geomean:',mu,' geosd:',gsd)
	rows = cur.execute("SELECT COUNT(*) from cudd_632ed69plus dmc WHERE (dmc.type='ps' OR dmc.type='bayes') AND (dmc.total_t < 1000)").fetchall()
	print('Total instances solved by top down',rows)
	rows = cur.execute("SELECT COUNT(*) from budmc dmc WHERE (dmc.type='ps' OR dmc.type='bayes') AND (dmc.total_t < 1000)").fetchall()
	print('Total instances solved by bottomup',rows)
	rows = cur.execute("SELECT "+dmc2TimeStr+"/"+dmc1TimeStr+",dmc2.declaredNodeCount from cudd_632ed69plus dmc1, budmc dmc2  WHERE "+titleStr+" AND (dmc1.name = dmc2.name) AND (dmc1.type = dmc2.type) AND ("+dmc1TimeStr+" < 1000) AND ("+dmc2TimeStr+" < 1000)").fetchall()
	from matplotlib import rc
	rc('text', usetex=True)
	colors=['green','blue','red','brown','orange']
	mkrs= ['o','^','+','x','D']
	speedups, nodecounts = zip(*rows)
	pl.scatter(nodecounts,speedups, color=colors[0], marker='x')
	pl.plot(np.unique(nodecounts), np.poly1d(np.polyfit(nodecounts,speedups, 1))(np.unique(nodecounts)), label='Best-Fit Line')
	#pl.errorbar([g[0] for g in gms],[g[1] for g in gms],[[g[1]-g[3] for g in gms],[g[2]-g[1] for g in gms]])
	# for g in gms:
		# pl.errorbar(g[0],g[2]-g[3],bottom=g[3])
	pl.rcParams['pdf.fonttype'] = 42
	pl.rcParams['ps.fonttype'] = 42	
	pl.legend(loc='upper right',fontsize=12)
	pl.tick_params(labelsize=16)
	pl.ylim(-5,280)
	pl.xlim(0,13000)
	pl.ticklabel_format(style='sci', axis='x', scilimits=(-3,4))
	pl.xlabel('\# nodes in Project-Join Tree', fontsize=20)
	pl.ylabel('Speedup', fontsize=20)
	#pl.title('Speedup of top-down sampling relative to bottom-up', fontsize=24)
	pl.axhline(y=1, color='r', linestyle='--')
	pl.savefig('result_analysis/graphs/dmc_632ed69+/sampling_comparison2.png',bbox_inches='tight')
	pl.show()
	
 
def calcSylvanCuddSpeedup(cur: Cursor, title:str):
	dmc1TimeStr = '(dmc1.totpreaddcomp_t + dmc1.addcomp_t)'
	dmc2TimeStr = '(dmc2.totpreaddcomp_t + dmc2.addcomp_t)'
	if title == 'all': #not including enc1.. do it separately if needed
		titleStr = '((dmc1.type=\'bayes\' and dmc2.type=\'bayes\') OR (dmc1.type=\'ps\' and dmc2.type=\'ps\') OR (dmc1.type=\'wb\' and dmc2.type=\'wb\'))'
	elif title == 'allunique':
		titleStr = '((dmc1.type=\'bayes\' and dmc2.type=\'bayes\') OR (dmc1.type=\'ps\' and dmc2.type=\'ps\') OR (dmc1.type=\'wb\' and dmc2.type=\'wb\'AND (dmc2.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')) ))'
	elif title=='bayes':
		titleStr = '(dmc1.type=\'bayes\' and dmc2.type=\'bayes\')'
	elif title == 'ps':
		titleStr = '(dmc1.type=\'ps\' and dmc2.type=\'ps\')'	
	elif title == 'wb':
		titleStr = '(dmc1.type=\'wb\' and dmc2.type=\'wb\')'
	elif title == 'db':
		titleStr = '((dmc1.type=\'bayes\' and dmc2.type=\'bayes\') OR (dmc1.type=\'ps\' and dmc2.type=\'ps\'))'
	rows = cur.execute("SELECT "+dmc1TimeStr+"/"+dmc2TimeStr+" from cudd_632ed69plus dmc1, sylvan_632ed69plus dmc2  WHERE "+titleStr+" AND (dmc1.name = dmc2.name) AND (dmc1.type = dmc2.type) AND ("+dmc1TimeStr+" < 1000) AND ("+dmc2TimeStr+" < 1000)").fetchall()
	print (geomean(next(zip(*rows))))

def failureAnalysis(cur):
	rows = cur.execute("SELECT COUNT(*) FROM cudd_632ed69plus dmc WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\'))))").fetchall()
	print('Total unique benchmarks:',rows)
	rows = cur.execute("SELECT COUNT(*) FROM cudd_632ed69plus dmc WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (dmc.total_t IS NULL)").fetchall()
	print('Unique Benchmarks where dmc (cudd) failed to generate 5000 samples',rows)
	rows = cur.execute("SELECT COUNT(*) FROM cudd_632ed69plus dmc WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (dmc.total_t IS NULL) AND (dmc.planner_t IS NULL)").fetchall()
	print('Unique Benchmarks where dmc (cudd) failed to generate 5000 samples and (because) no pj-tree was computed',rows)
	rows = cur.execute("SELECT COUNT(*) FROM cudd_632ed69plus dmc, waps WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (waps.type=dmc.type) AND (waps.name=dmc.name) AND (waps.total_t <1000) AND (dmc.total_t IS NULL)").fetchall()
	print('Total unique benchmarks where Waps succeeded but dmc failed',rows)
	rows = cur.execute("SELECT COUNT(*) FROM cudd_632ed69plus dmc, waps WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (waps.type=dmc.type) AND (waps.name=dmc.name) AND (waps.total_t <1000) AND (dmc.total_t IS NULL) AND (dmc.planner_t IS NULL)").fetchall()
	print('Total unique benchmarks where Waps succeeded but dmc failed and (because) pj-tree was not constructed',rows)
	rows = cur.execute("SELECT COUNT(*) FROM cudd_632ed69plus dmc WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (dmc.planner_t < 1000)").fetchall()
	print('Total unique benchmarks where pj-tree was constructed',rows)
	denom = rows[0][0]
	rows = cur.execute("SELECT COUNT(*)/"+str(denom)+".0 FROM cudd_632ed69plus dmc WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (dmc.planner_t < 1000) AND (dmc.addcomp_t < 1000)").fetchall()
	print('fraction of unique pj-tree constructed where add compilation succeeded',rows)
	rows = cur.execute("SELECT COUNT(*), COUNT(*)/"+str(denom)+".0 FROM cudd_632ed69plus dmc WHERE (dmc.type='bayes' OR dmc.type='ps' OR (dmc.type='wb'AND (dmc.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (dmc.planner_t < 1000) AND (dmc.addcomp_t < 1000) AND (dmc.total_t < 1000)").fetchall()
	print('fraction of unique pj-tree constructed which were fully solved',rows)
	
def main():
	con, cur = initDB()
	#not including enc1 in all. Do it separately if needed.
	#rows = cur.execute("SELECT COUNT(*) from cudd_632ed69plus dmc, waps WHERE (dmc.name=waps.name) AND (dmc.type=waps.type) AND (waps.type='bayes' OR waps.type='ps' OR (waps.type=\'wb\' AND (waps.name NOT IN (SELECT waps2.name FROM waps waps2 WHERE waps2.type=\'bayes\' OR waps2.type=\'ps\')))) AND (dmc.total_t < 1000 AND waps.total_t < 1000)").fetchall()
	#print(rows)
	#sys.exit(1)
	print("Geometric Mean of Speedups relative to WAPS:")
	for t in ['totaltime', 'compiletime']:
		lib='cudd'
		for title in ['db','wb','all','allunique']:
			print("	"+t+" "+lib+" "+title+":",end=" ")
			calcSpeedups(cur,t=='totaltime',lib,title)
	print("Number of benchmarks finished:")
	for t in ['totaltime', 'compiletime']:
		lib='cudd'
		for title in ['db','wb','all','allunique']:
			print("	"+t+" "+lib+" "+title+":",end=" ")
			calcFinished(cur,t=='totaltime',lib,title)
	print("Number of benchmarks fastest on (total time):")
	lib='cudd'
	for title in ['db','wb','all','allunique']:
		print("	"+lib+" "+title+":",end=" ")
		calcFastest(cur,lib,title)
	
	for title in [['db','db',1945],['wb','wb', 773]]:
		plotCactus(cur,title[0],title[1],title[2])
	# plotSamplingComparison()
	
	print("Geometric Mean of Compile Time Sylvan Speedups relative to CUDD:")
	for title in ('db','wb','all','allunique'):
		print(title,end=': ')
		calcSylvanCuddSpeedup(cur, title)
	
	plotPar2TWs(cur,'allunique', 10,np.mean)
	
	print('Top vs bottom')
	plotSamplingComparison2(cur,'all')
	
	failureAnalysis(cur)
	#plotPar2TWs(cur,'bayes', 10,np.mean)
	#plotPar2TWs(cur,'ps', 10,np.mean)
	#plotPar2TWs(cur,'wb', 10,np.mean)
	
	# for t in ['totaltime', 'compiletime']:
	# 	for lib in ['cudd', 'sylvan']:
	# 		for title in ['bayes','ps','enc1','all']:
	# 			plotSpeedupTWs(cur,t=='totaltime',lib,title)
	
	#print('DPS(cudd) vs WAPS')
	#calcSpeedups(cur, False, 'cudd','wb')
	#calcFinished(cur, False, 'cudd','wb')
	#calcFastest(cur, 'cudd','wb')
	#print('CUDD vs Sylvan')
	#calcSylvanCuddSpeedup(cur, 'wb')
	
	closeAll(con, cur)

if __name__=="__main__":
   main()
