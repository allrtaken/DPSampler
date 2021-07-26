#ifndef ADDSAMPLER_H_
#define ADDSAMPLER_H_

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "RandomBits.hpp"
//#include "sampler/SamplerNodeFactory.hpp"
#include <gmpxx.h>

#include "../logic.hh"
#include "../Dd.hh"

using std::vector;
using std::unordered_set;
using std::string;
using std::unordered_map;
//using Sampler::SamplerNodeFactory;

namespace Sampler{
#if SAMPLE_NUM_TYPE == 2
	typedef mpf_class WtType;
#else
	typedef Float WtType;
#endif

class Asmt{
	public:
		Asmt(Int nVars_): bits(nVars_,-1){  //-1 is unset default value
		}
		void setNoCheck(bool bit, Int i);
		void set(bool bit, Int i);
		bool get(Int i);
		void clear();
		void printAsmt();
		bool isSet(Int i);
	private:
		vector<int8_t> bits; //not bool becz need to have 3rd 'unset' value (defined to be -1)
};

/*
sign of cmprsdLvl stores whether nodevar is part of sample-set (+ve sign) or if it is previously assigned (-ve sign)
if its in sample-set (+ve), we store cmprsdLvl in var cmprsdLvl, other wise we store cnfvarindex.
We use least significant bit of cmprsdLvl to store whether the SamplerNode has been visited previously (by inplacecofactor).
inplacecofactor sets this bit after visiting the node, and inplaceuncofactor resets it
In order to accommodate this bit, we left shift (multiply by two) cmprsdLvl/cnfvarindex before storing in cmprsdLvl
*/
class SamplerNode{
	public:
	SamplerNode(WtType, SamplerNode*, SamplerNode*, Int);
	WtType wt;
	SamplerNode* t;
	SamplerNode* e;
	Int cnfVarID;
	//if node var is part of sampleset, auxvar stores wt*factor of 'then' branch, 
	//if node var is previously assigned, auxvar stores inherited cmprsdLvl from appropriate child node
	WtType auxVar; 
	// DdNode* dnode;
};

class ADDSampler{
	public:
		ADDSampler(const JoinNode* root_, const Cudd* mgr_, Int nTotalVars_, Int nApparentVars, unordered_map<Int,Int> c2DVarMap, vector<Int> d2CVarMap, 
			unordered_map<Int, Number> litWeights_, Set<Int> freeVars, bool checkAsmts_);
		
		void buildDataStructures();
		Asmt& drawSample();
		void writeAsmtToFile(FILE* ofp);
		void sampleAndWriteToFile(FILE* ofp, Int nSamples, bool printProgress);
		//Int getNumVars();
		
	private:
		void createAuxStructures(const JoinNode*);
		void createAuxStructure(const JoinNode*);
		void createSamplingDAGs(const JoinNode*);
		SamplerNode* createSamplingDAG(DdNode* node, Dd* dd, unordered_map<DdNode*, SamplerNode*>&);
		WtType sampleFromADD(SamplerNode*, vector<Int>&, Set<Int>&);
		void drawSample_rec(const JoinNode*);		
		double getAsmtVal(Dd*);
		//WtType inplaceCofactor(SamplerNode*, vector<Int>&);
		Int inplaceCofactor(SamplerNode* sNode, Int* sampleCNFVarIDs, Int** currVarPtr);
		void inplaceUnCofactor(SamplerNode*);
		//WtType computeFactor(Int startCmprsdLvl, Int endCmprsdLvl, bool tE, vector<Int>&);
		WtType computeFactor(Int* sampleVarCNFIDs, Int currCNFVarID, Int** currVarPtr, bool tE);
		void sampleMissing(Set<Int>& notSampled);
		void seekForward(Int* sampleVarCNFIDs, Int currCNFVarID, Int** currVarPtr);

		bool checkAsmts;

		Int nTotalVars, nApparentVars;
		Asmt* t;
		DdNode** sVars;
		RandomBits rb;
		//SamplerNodeFactory *snFactory;

		const JoinNode* jtRoot;
		unordered_map<Int, Int> cnfVarToDdVarMap;
		vector<Int> ddVarToCnfVarMap;
		ADD assignedVarsCube;
		int numAssigned;
		Cudd mgr;
		//unordered_map<Int, Float> litWeights;
		vector<WtType> litWts;
		unordered_map<const JoinNode*, SamplerNode*> rootSNMap;
		unordered_map<DdNode*, SamplerNode*> nodeMap;
		Set<Int> freeVars;
		//Float testNum, testDen;

		//TimePoint startTime, allAuxStructsTime, allDAGsTime, allSampleGenTime;
};
}//end namespace
#endif