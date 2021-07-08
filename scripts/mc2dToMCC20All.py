from runOnAll import runOnAll
from sys import argv

inDir = argv[1]
outDir = argv[2]

cmdList = ['python3.7','mc2dToMCC20.py',]
runOnAll(inDir=inDir,inExt='cnf',cmdList=cmdList,inFIndex=2,outDir=outDir,outExt='cnf',outFIndex=3,cmdOutFile='convert.log',skip=True)