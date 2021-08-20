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

#include "../../libraries/sylvan/src/sylvan_mtbdd.h"
#include "../../libraries/sylvan/src/sylvan_int.h"

using std::vector;
using std::unordered_set;
using std::string;
using std::unordered_map;
//using Sampler::SamplerNodeFactory;

using sylvan::MTBDD;
using sylvan::mtbdd_isleaf;
using sylvan::mtbdd_gethigh;
using sylvan::mtbdd_getlow;
using sylvan::mtbdd_getvar;
using sylvan::MTBDD_GETNODE;
using sylvan::mtbddnode_t;

namespace Sampler{
#if SAMPLE_NUM_TYPE == 2
	typedef mpf_class WtType;
#else
	typedef Float WtType;
#endif
extern bool usingCUDD;

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

class SamplerUtils{
	public:
	double_t getLog(mpz_t inputNum);
	void fp_epsilon (mpf_t eps, uint32_t prec);
	void fp_log_m1 (mpf_t lg, const mpf_t z, uint32_t prec);
	void fp_log2 (mpf_t l2, uint32_t prec);
	void fp_log (mpf_t lg, const mpf_t z, uint32_t prec);

};

//sign of cnfvarid denotes if the node has been visited previously or not
//auxvar technically not required if we normalize wts to 1. But currently memory is not an issue and
//storing borrowed cnfvarid (see below) is helpful
class SamplerNode{
	public:
	SamplerNode(WtType, SamplerNode*, SamplerNode*, Int);
	WtType wt;
	SamplerNode* t;
	SamplerNode* e;
	Int cnfVarID;
	//auxvar stores wt of then branch for fast random sampling. In case of non-samplevarnodes it stores borrowed
	//cnfvarid from child
	WtType auxVar; 
	// DdNode* dnode;
};

class ADDSampler{
	public:
		ADDSampler(const JoinNode* root_, const Cudd* mgr_, Int nTotalVars_, Set<Int> apparentVars_, Int ddNodeCount_, unordered_map<Int,Int> c2DVarMap, vector<Int> d2CVarMap, 
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
		SamplerNode* createSamplingDAG(DdNode* node, Dd* dd);
		SamplerNode* createSamplingDAG(MTBDD node, Dd* a);
		WtType sampleFromADD(SamplerNode*, Set<Int>&);
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
		Int ddNodeCount;
		Asmt* t;
		DdNode** sVars;
		RandomBits rb;
		//SamplerNodeFactory *snFactory;

		const JoinNode* jtRoot;
		unordered_map<Int, Int> cnfVarToDdVarMap;
		vector<Int> ddVarToCnfVarMap;
		// ADD assignedVarsCube;
		int numAssigned;
		Cudd mgr;
		//unordered_map<Int, Float> litWeights;
		vector<WtType> litWts;
		unordered_map<const JoinNode*, SamplerNode*> rootSNMap;
		unordered_map<DdNode*, SamplerNode*> nodeMap_cudd;
		unordered_map<mtbddnode_t, SamplerNode*> nodeMap_sylvan;
		Set<Int> freeVars;
		Set<Int> apparentVars;
		Set<Int> allProjVars;
		//Float testNum, testDen;

		//TimePoint startTime, allAuxStructsTime, allDAGsTime, allSampleGenTime;
};
}//end namespace
#endif