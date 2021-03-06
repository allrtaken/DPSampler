import sys,os

cnf = sys.argv[1]
cs = sys.argv[2]
ns = 1000

assert(os.environ.has_key('DPSAMPLER'))

if len(sys.argv) > 3:
	ns = int(sys.argv[3])
w = 0.5
if len(sys.argv) > 4:
	w_check = float(sys.argv[4])
	if w_check > 0:
		w = w_check
sf = 'samples.txt'
if len(sys.argv) > 5:
	sf = sys.argv[5]

'''
cnf2 = '/tmp/wtd'+'_wt='+str(w)+'.cnf'
f1 = open(cnf,'r')
f2 = open(cnf2,'w')
addWts = True
for line in f1:
	f2.write(line)
	if line.startswith('p cnf'):
		nVars = int(line.split()[2])
		print 'nVars=',nVars
		wts = ''
		for i in range(2*nVars):
			wts += str(w)+' '
	if line.startswith('c weights'):
		addWts = False
if addWts:
	f2.write('c weights '+wts+'\n')
f1.close()
f2.close()
cnf = cnf2
'''

cmd1 = os.environ['DPSAMPLER']+'/lg/build/lg "lg/solvers/flow-cutter-pace17/flow_cutter_pace17 -s 1234567 -p 100" < '+cnf+' > tree.tmp'
cmd2 = os.environ['DPSAMPLER']+' DMC/dmc --cf='+cnf+' --cs='+cs+' --sf='+sf+'< tree.tmp'
cmd3 = os.environ['DPSAMPLER']+'lg/build/lg "'+os.environ['DPSAMPLER']+'lg/solvers/flow-cutter-pace17/flow_cutter_pace17 -s 1234567 -p 100" < '+cnf+' | '+os.environ['DPSAMPLER']+'dmc/dmc --pw=0 --cf='+cnf+' --lc=1 --wc=1 --cs='+cs+' --sf='+sf+' --ns='+str(ns)
print cmd3
os.system(cmd3)
#print cmd2
#os.system(cmd2)
'''
$DPSAMPLER/lg/build/lg "$DPSAMPLER/lg/solvers/flow-cutter-pace17/flow_cutter_pace17 -s 1234567 -p 100" < test.cnf | $DPSAMPLER/dmc/dmc --pw 0 --cf test.cnf --cs c
$DPSAMPLER/lg/build/lg "$DPSAMPLER/lg/solvers/flow-cutter-pace17/flow_cutter_pace17 -s 1234567 -p 100" < converted_benchmarks/enc1/avgdeg_3_016_1.txt.cnf | $DPSAMPLER/dmc/dmc --pw 0 --cf converted_benchmarks/enc1/avgdeg_3_016_1.txt.cnf --cs c --jp a --dp s --lc 0 --tc 0 
'''
