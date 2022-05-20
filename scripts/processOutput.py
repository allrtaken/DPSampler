from sqlite3.dbapi2 import Cursor
import sys, os
from outputFileProcessors import processJTFile, processDPSFile, processWAPSFile, processTimeFile, processBUDPSFile
import sqlite3 as sq
from contextlib import closing
from benchmark_names import bayes_names,ps_names,enc1_names, wb_names

# I manually worked with following functions in the python interpreter to see that everything was working as expected
# therefore no 'main' function here.
# some commands:
# from processOutput import initDB, getTableList, getTableSchema, closeAll
# conn, cursor = initDB()
# getTableList(cursor)
# getTableSchema(cursor, "waps")
# conn.commit()
# closeAll(conn, cursor)
# cursor.execute("ALTER TABLE sylvan_632ed69plus ADD mem INTEGER")
# cursor.execute("CREATE TABLE cudd_632ed69plus (name TEXT, type TEXT, jt_width INTEGER, planner_t REAL, totpreaddcomp_t REAL, addcomp_t REAL, sampcomp_t REAL, sampgen_t REAL, total_t REAL, PRIMARY KEY(name, type))")
# cursor.execute("CREATE TABLE budmc (name TEXT, type TEXT, jt_width INTEGER, declaredNodeCount INTEGER, sampgen_t REAL, total_t REAL, mem INTEGER, PRIMARY KEY (name, type))")

def initDB():
	conn = sq.connect('dpsampling.db')
	print("DB connection total changes: ",conn.total_changes)
	cursor = conn.cursor()
	return conn,cursor

def closeAll(cn, cr):
	cr.close()
	cn.close()

def getTableList(cursor):
	cursor.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;")
	tabs = (cursor.fetchall())
	print(tabs)

def getTableSchema(cursor, tab):
	s = "PRAGMA table_info('"+tab+"')"
	# print(s)
	for row in cursor.execute(s).fetchall():
		print(row)

def fillBenchNames(cursor: Cursor):
	flist = None
	# flists = [[bayes_names,'bayes'],[ps_names,'ps'],[enc1_names,'enc1']]
	flists = [[wb_names,'wb']]
	for flist,title in flists:
		tab_name = 'filelist'
		print("Inserting filenames and indices of "+title+" benchmarks into table "+tab_name+"..")
		skipped = []
		for i,name in enumerate(flist):
			rows = cursor.execute("SELECT * FROM "+tab_name+" WHERE name=? AND type=?",(name,title)).fetchall()
			assert len(rows) <= 1 #name,type is primary key
			if len(rows) == 1:
				if rows[0][1] != i:
					print("ERROR: Table "+tab_name+" already contains a record with name=",rows[0][0]," and type=",rows[0][1]," and index=",rows[0][2])
					print(".. when trying to insert new record with name=",name," and index=",i,"with type",title)
					print("Stopping insertions..")
					break
				else:
					skipped.append(i)
					continue
			cursor.execute("INSERT INTO "+tab_name+" VALUES (?,?,?)",(name,title,i))
		print("Skipped following",len(skipped),"indices because records were already present:",skipped)
		rows = cursor.execute("SELECT COUNT(1) FROM "+tab_name+"").fetchall()
		print("Finished inserting file names in "+tab_name+". #records in "+tab_name+" is ", rows[0][0])

def fillDMC(cursor: Cursor, dDir,lib):
	params = [[1080,'bayes','bayes'],[896,'pseudoweighted','ps'],[270,'enc1','enc1']]
	params = [[773,'waps_bench','wb']]
	for param in params:
		print("Doing "+param[1])
		for i in range(param[0]):
			fp = open(dDir+'/'+lib+'/output/'+param[1]+'/dps_array_'+str(i)+'.out','r')
			fName, joinTreeWidth, plannerTime, preADDCompTime, ADDCompTime, SampCompTime, SampGenTime, totTime = processDPSFile(fp)
			fp.close()
			fp = open(dDir+'/'+lib+'/timing/'+param[1]+'/dps_array_'+str(i)+'.time','r')
			userT, sysT, elapsedT, mem = processTimeFile(fp)
			fp.close()
			tup = (fName,param[2],joinTreeWidth, plannerTime, preADDCompTime, ADDCompTime, SampCompTime, SampGenTime, totTime, mem)
			cursor.execute("INSERT INTO "+lib+"_632ed69plus VALUES (?,?,?,?,?,?,?,?,?,?)",tup)
		rows = cursor.execute("SELECT COUNT(1) FROM "+lib+"_632ed69plus").fetchall()
		print("Finished inserting file names in "+lib+"_632ed69plus. #records in "+lib+"_632ed69plus is ", rows[0][0])

def fillBUDMC(cursor: Cursor, dDir):
	params = [[1080,'bayes','bayes'],[896,'pseudoweighted','ps']]
	for param in params:
		print("Doing "+param[1])
		for i in range(param[0]):
			fp = open(dDir+'/output/'+param[1]+'/dps_array_'+str(i)+'.out','r')
			fName, joinTreeWidth, declaredNodeCount, SampGenTime, totTime = processBUDPSFile(fp)
			fp.close()
			fp = open(dDir+'/timing/'+param[1]+'/dps_array_'+str(i)+'.time','r')
			userT, sysT, elapsedT, mem = processTimeFile(fp)
			fp.close()
			tup = (fName,param[2],joinTreeWidth, declaredNodeCount, SampGenTime, totTime, mem)
			cursor.execute("INSERT INTO budmc VALUES (?,?,?,?,?,?,?)",tup)
		rows = cursor.execute("SELECT COUNT(1) FROM budmc").fetchall()
		print("Finished inserting file names in budmc. #records in budmc is ", rows[0][0])

def fillWAPS(cursor: Cursor, wDir):
	params = [[1080,'bayes','bayes'],[896,'pseudoweighted','ps'],[270,'enc1','enc1']]
	params = [[773,'waps_bench','wb']]
	for param in params:
		print("Doing "+param[1])
		for i in range(param[0]):
			fp = open(wDir+'/old_3_output/'+param[1]+'/waps_array_'+str(i)+'.out','r')
			fName, compileTime, totTime = processWAPSFile(fp)
			fp.close()
			fp = open(wDir+'/old_3_timing/'+param[1]+'/waps_array_'+str(i)+'.time','r')
			userT, sysT, elapsedT, mem = processTimeFile(fp)
			fp.close()
			tup = (fName,param[2], compileTime, totTime, mem)
			cursor.execute("INSERT INTO waps VALUES (?,?,?,?,?)",tup)
		rows = cursor.execute("SELECT COUNT(1) FROM waps").fetchall()
		print("Finished inserting file names in waps. #records in waps is ", rows[0][0])
	
# def main():
# 	initDB()
# 	#populateDatabase()

# if __name__=="__main__":
#    main()
