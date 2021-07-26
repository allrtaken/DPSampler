#include <sys/types.h>
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <random>
#include <cassert>

#include "ADDSampler.hpp"

using std::cout;
using std::endl;
using std::sort;
using Sampler::Asmt;
using Sampler::ADDSampler;
using Sampler::SamplerNode;
using Sampler::WtType;

extern Int samplePrec;
extern const Float INF;

void Asmt::set(bool bit, Int i){
	assert(bits.at(i)==-1); //sampling is generally not the bottleneck so leaving asserts in place
	bits[i] = bit;
}
bool Asmt::isSet(Int i){
	return (bits.at(i)!=-1);
}
void Asmt::setNoCheck(bool bit, Int i){
	printf("\nSet called on state %p at pos %lld!\n",this,i);
	bits[i] = bit;
}

void Asmt::printAsmt(){
	cout<<"Asmt:\n";
	for(Int i = 0; i<bits.size();i++){
		printf("%d",bits[i]);
	}
	cout<<"\n";
}

bool Asmt::get(Int i){
	assert(bits.at(i)!=-1); //sampling is generally not the bottleneck so leaving asserts in place
	return bits[i]==1;
}
void Asmt::clear(){
	for (Int i = 0; i<bits.size();i++){
		bits[i] = -1;
	}
}

Sampler::SamplerNode::SamplerNode(WtType wt_, SamplerNode* t_, SamplerNode* e_, Int cnfVarID_): wt(wt_), t(t_), e(e_), 
	cnfVarID(cnfVarID_), auxVar(1) {}

ADDSampler::ADDSampler(const JoinNode *root_, const Cudd* mgr_ , Int nTotalVars_, Int nApparentVars_, unordered_map<Int,Int> c2DVarMap, vector<Int> d2CVarMap,
	unordered_map<Int, Number> litWeights_, Set<Int> freeVars_,	bool checkAsmts_): cnfVarToDdVarMap(c2DVarMap), ddVarToCnfVarMap(d2CVarMap), 
	litWts(litWeights_.size()*3/2+3), freeVars(freeVars_), checkAsmts(checkAsmts_), jtRoot(root_), mgr(*mgr_), nTotalVars(nTotalVars_),
	nApparentVars(nApparentVars_){
	
	rb.SeedEngine();
	rb.SeedEngine2();
	litWts.at(0) = -1;
	litWts.at(1) = -1;
	litWts.at(2) = -1;
	
	#ifndef SAMPLE_NUM_TYPE
		#define SAMPLE_NUM_TYPE 1
		util::printComment("Macro SAMPLE_NUM_TYPE not defined. Setting to log-counting (1)..");
	#endif 

	#if SAMPLE_NUM_TYPE == 0
		util::printComment("Computing sampling weights using Floats (See SAMPLE_NUM_TYPE compile-time option).");
	#elif SAMPLE_NUM_TYPE == 1
		util::printComment("Computing sampling weights using log-Counting (See SAMPLE_NUM_TYPE compile-time option).");
	#elif SAMPLE_NUM_TYPE == 2
		util::printComment("Computing sampling weights using MPF(GMP) with precision set to "+to_string(samplePrec)+" bits (See SAMPLE_NUM_TYPE compile-time option).");
		mpf_set_default_prec(samplePrec);
	#else
		util::printComment("Macro SAMPLE_NUM_TYPE not recognized. Setting to default (0)..");
		#define SAMPLE_NUM_TYPE 0
		util::printComment("Computing sampling weights using Floats (See SAMPLE_NUM_TYPE compile-time option).");	
	#endif
	util::printComment("Note: This determines how the weights for the sampler are computed. "
	"This is independent of logCounting or multiprecision flags of dmc which are only used for ADD-construction.");
	for(Int i = 1; i<=nTotalVars;i++){
		
		#if SAMPLE_NUM_TYPE == 0
			litWts.at(3*i) = litWeights_.at(i).fraction;
			litWts.at(3*i+1) = litWeights_.at(i).fraction+litWeights_.at(-i).fraction;
			litWts.at(3*i+2) = litWeights_.at(-i).fraction;
		#elif SAMPLE_NUM_TYPE == 1
			litWts.at(3*i) = log10l(litWeights_.at(i).fraction);
			litWts.at(3*i+1) = log10l(litWeights_.at(i).fraction+litWeights_.at(-i).fraction);
			litWts.at(3*i+2) = log10l(litWeights_.at(-i).fraction);
		#elif SAMPLE_NUM_TYPE == 2
			/*Cast to double because MPF does not have constructors for long double*/
			litWts.at(3*i) = (double)litWeights_.at(i).fraction;
			litWts.at(3*i+1) = (double) (litWeights_.at(i)+litWeights_.at(-i)).fraction;
			litWts.at(3*i+2) = (double) (litWeights_.at(-i)).fraction;
			//litWts.at(2*i).canonicalize();
		#endif
	}
	// cout<<"\n";
	for (Int i = 0; i<nApparentVars; i++){
		assert(mgr.ReadPerm(i)==i);
		//cout<<i<<" "<<mgr.ReadPerm(i)<<"\n";
		//this is because dpmc maintains mapping between cnfvars and ddvars instead of using cudd_addpermute
		//so index and levels of all ddvars should be same. 
		//this is different from tracesampler which used permute
		//therefore all code here is written under assumption that index and levels are same for all ddvars
	}
	assert(nApparentVars+freeVars.size()==nTotalVars);
	numAssigned = 0;
}

void ADDSampler::buildDataStructures(){
	t = new Asmt(nTotalVars+1); //cnfvarids are indexed from 1.
	//Cudd_AutodynDisable(a.dd);
	//printComment("#nodes in JoinTree:"+to_string(jtRoot->getNodeCount()));
	util::printComment("Building aux structures..  ",0,0);
	createAuxStructures(jtRoot);
	// allAuxStructsTime = util::getTimePoint();
	util::printComment("Built all aux structures!",0,1,false);
	util::printComment("Building sampling DAGs..  ",0,0);
	createSamplingDAGs(jtRoot);
	util::printComment("Built sampling DAGs!",0,1,false);
	util::printComment("Finished building all Datastructures!");
}

void ADDSampler::createAuxStructure(const JoinNode* jNode){
	//compressedPerm only includes levels for variables that are to be sampled, i.e. the vars in projectedVars of joinnode
	//that means that for root add, it will include all vars in the add
	class Indexer{
		public:
		Int level;
		Int index;
		Indexer(Int index_, Int level_): level(level_), index(index_){}
		bool operator < (const Indexer& str) const{
        	return (level < str.level);
    	}
	};
	assert(jNode!=NULL);
	Set<Int> pVars = jNode->projectionVars;
	Int numPVars = pVars.size();
	if(numPVars == 0){
		return;
	}
	Dd* a = (Dd*)jNode->dd;
	assert(a!=NULL);
	if (jNode==jtRoot){
		a->precompileSampleDAG = true;  // no assignments at root so can precompile
	} else{
		a->precompileSampleDAG = true; //precompile all for new top-down sampling method
		// for older bottom up sampling cant precompile unless restrictions on variable ordering.
		//so had set all non-root nodes to false. Now can precompile everything
	}
	vector<Indexer> ind_level_map;
	//generate perm o(n)
	for(auto p: pVars){
		ind_level_map.emplace_back(p,cnfVarToDdVarMap.at(p));
	}
	a->sampleStructDepth = numPVars + 1; //since sample structure should include leaves
	
	// cout<<"Sorting..\n";
	//sort nlogn
	sort(ind_level_map.begin(),ind_level_map.end());
	//compress and insert o(n) 
	
	a->compressedInvPerm.resize(numPVars*2);
	a->sampleVarCNFIDs = new Int[numPVars+2]; //+2 sentinels
	a->sampleVarCNFIDs[0] = MIN_INT;
	for(int i = 0; i< numPVars; i++){
		a->sampleVarCNFIDs[i+1] = ind_level_map.at(i).index;
		a->compressedPerm[ind_level_map.at(i).level] = i; //cmprsdperm stores mapping from ddvarid to cmprsdlevl
		a->compressedInvPerm.at(2*i) = ind_level_map.at(i).index;
		a->compressedInvPerm.at(2*i+1) = ind_level_map.at(i).level;
	}
	a->sampleVarCNFIDs[numPVars+1] = MAX_INT;
	a->curVarPtr = a->sampleVarCNFIDs+1;
	//cmprsdperm from ddvarindex to cmprsdlvl
	//cmprsdinvprm from cmprsdlvl to cnfvarindex,ddvarindx
}

void ADDSampler::createAuxStructures(const JoinNode* jNode){
	createAuxStructure(jNode);
	for (auto child: jNode->children){
		createAuxStructures(child);
	}
}

void ADDSampler::createSamplingDAGs(const JoinNode* jNode){
	//cout<<"Creating sampling structure..\n";
	assert(jNode!=NULL);
	if (jNode->projectionVars.size()==0){
		//cout<<"Sampling DAG skipped because projvar size is 0\n";
		return;
	}
	Dd* a = (Dd*)jNode->dd;
	assert(a!=NULL);
	if(!a->precompileSampleDAG) {
		//cout<<"Not precompiling\n";
	} else{
		//unordered_map<DdNode*, SamplerNode*> nodeMap;
		SamplerNode* sNode = createSamplingDAG(a->cuadd.getNode() ,a, nodeMap);
		rootSNMap[jNode] = sNode;
	}
	for (auto child: jNode->children){
		createSamplingDAGs(child);
	}
}

SamplerNode* ADDSampler::createSamplingDAG(DdNode* node, Dd* a, unordered_map<DdNode*, SamplerNode*>& nodeMap){
	assert(node!=NULL);
	if (nodeMap.find(node)!=nodeMap.end()){ 
		SamplerNode* snode = nodeMap.at(node);
		return snode;
	}
	SamplerNode *sNode = NULL, *sNode_t = NULL, *sNode_e = NULL;
	WtType wt = 0;
	Int cnfVarID;
	if (!Cudd_IsConstant(node)){
		Int nodeInd = node->index;
		cnfVarID = ddVarToCnfVarMap.at(nodeInd);
		sNode_t = createSamplingDAG(Cudd_T(node), a, nodeMap);
		sNode_e = createSamplingDAG(Cudd_E(node), a, nodeMap);
	} else{
		if (logCounting){
			#if SAMPLE_NUM_TYPE == 0 || SAMPLE_NUM_TYPE == 2
				wt = (double) exp10l(Cudd_V(node));
			#else
				wt = Cudd_V(node); 
			#endif
		} else {
			#if SAMPLE_NUM_TYPE == 0 || SAMPLE_NUM_TYPE == 2
				wt = Cudd_V(node);
			#else
				wt = log10l(Cudd_V(node));
				// if (Cudd_V(node)==0){
				// 	cout<<"log of 0 is "<<wt<<" and the constant is "<<-INF<<" and equality is "<<(wt==-INF)<<"\n";
				// }
			#endif
		}
		
		cnfVarID = MAX_INT; // leaf level should be equal to size of projectable vars, since cmprsdlvls are 0 indexed
	}
	sNode = new SamplerNode(wt,sNode_t,sNode_e,cnfVarID);
	// sNode->dnode = node;
	nodeMap.emplace(node, sNode);
	return sNode;
}

WtType ADDSampler::computeFactor(Int* sampleVarCNFIDs, Int currCNFVarID, Int** currVarPtr, bool tE){
	// cout<<"currCNFVARID in computeFactor:"<<currCNFVarID<<" "<<std::flush;
	WtType factor = tE? litWts.at(currCNFVarID*3):litWts.at(currCNFVarID*3+2);
	(*currVarPtr)--; // see semantics (invaraints) of inpcf
	while (**currVarPtr != currCNFVarID){
		if (**currVarPtr == MIN_INT){
			cout<<"ERROR: current cnfvarid "<<currCNFVarID<<" not found while computing factor. Exiting..\n";
			exit(1);
		}
		#if SAMPLE_NUM_TYPE == 0 || SAMPLE_NUM_TYPE == 2			
			factor *= litWts.at(3*(**currVarPtr)+1);
		#else
			factor += litWts.at(3*(**currVarPtr)+1);
		#endif
		(*currVarPtr)--;
	}
	#if SAMPLE_NUM_TYPE == 0 || SAMPLE_NUM_TYPE == 2
		assert(factor!=0);
	#else
		assert(factor!=-INF);
	#endif
	return factor;
}

void ADDSampler::seekForward(Int* sampleVarCNFIDs, Int currCNFVarID, Int** currVarPtr){
	while (**currVarPtr != currCNFVarID && **currVarPtr != MAX_INT){
		(*currVarPtr)++;
	}
}

//TODO: seek forward to topmost samplevar or maxint
//computefactor should decrement first then mulitply
// nonsamplevar nodes should store original cnfvarid in auxvar and borrowed id in normal
// uncofactor should restore correctcnfvarids for non-sample-varnodes
// it needs to detect somehow that node is not-samplervar

//invariant: currvarptr points to topmost samplevarcnfid (or maxint) in the subtree in samplecnfvarids when the function returns
//invariant: for samplevarnodes, wt variable contains cumulative and scaled weight of subtree below that node
//invariant: for non-samplevarnodes, wt variable containts cumulative weight of subtree below topmost samplevar below node
// 			 in other words, for non-samplevarnodes, wt may not be scaled
Int ADDSampler::inplaceCofactor(SamplerNode* sNode, Int* sampleCNFVarIDs, Int** currVarPtr){
	// cout<<"Inside cofactor with value at currVarPtr at "<<**currVarPtr<<" .. "<<std::flush;
	if (sNode->t==NULL){
		//snode should be leaf
		assert(sNode->e==NULL);
		// assert(sNode->wt == Cudd_V(sNode->dnode));
		seekForward(sampleCNFVarIDs,MAX_INT,currVarPtr);
		return MAX_INT;
	}
	if (sNode->cnfVarID < 0){ // node has been visited previously
		if (!t->isSet(-(sNode->cnfVarID))){
			//if node is samplevar then use its cnfvarid
			seekForward(sampleCNFVarIDs,-(sNode->cnfVarID),currVarPtr);
			return -(sNode->cnfVarID);
		} else{
			//else get borrowed cnfvarid stored in auxvar
			Int borrowedCNFVARID;
			#if SAMPLE_NUM_TYPE != 2
				borrowedCNFVARID = sNode->auxVar;
			#else
				borrowedCNFVARID = sNode->auxVar.get_d(); //returns double which will be automatically casted to Int
			#endif
			seekForward(sampleCNFVarIDs,borrowedCNFVARID,currVarPtr);
			return borrowedCNFVARID;
		}
	}
	Int ret_CnfVarID;
	//if(cond is true), var is projectable at this node, therefore not previously assigned (will be assigned now)
	if (!t->isSet(sNode->cnfVarID)){
		// cout<<"Searching cnfVar "<<sNode->cnfVarID<<" "<<std::flush;
		seekForward(sampleCNFVarIDs,sNode->cnfVarID,currVarPtr);
		if (**currVarPtr == MAX_INT){
			cout<<"ERROR: current node's cnfvarid "<<sNode->cnfVarID<<" not found while doing inplacecofactor. Exiting..\n";
			exit(1);
		}
		// cout<<"Found cnfVar\n";
		Int childCnfVarID;
		#if SAMPLE_NUM_TYPE != 1
			childCnfVarID = inplaceCofactor(sNode->t, sampleCNFVarIDs, currVarPtr);
			WtType tFactor = computeFactor(sampleCNFVarIDs, sNode->cnfVarID, currVarPtr, true);
			WtType tWt = sNode->t->wt;
			if (tWt != 0) {
				sNode->auxVar = tWt * tFactor;
			} else {
				sNode->auxVar = 0;
			}
			childCnfVarID = inplaceCofactor(sNode->e, sampleCNFVarIDs, currVarPtr);
			WtType eFactor = computeFactor(sampleCNFVarIDs, sNode->cnfVarID, currVarPtr, false);
			WtType eWt = sNode->e->wt;
			if (eWt != 0) {
				sNode->wt = sNode->auxVar + eWt * eFactor;
			} else {
				sNode->wt = sNode->auxVar;
			}
		#else
			bool thenIsZero = false;
			// cout<<"Calling then-inpcf for curNodeID "<<sNode->cnfVarID<<".. "<<std::flush;
			childCnfVarID = inplaceCofactor(sNode->t, sampleCNFVarIDs, currVarPtr);
			//pulled this out of if-condition below, because computeFactor resets the currvarptr to cnfvarid
			//which is needed for else cond also.
			WtType tFactor = computeFactor(sampleCNFVarIDs, sNode->cnfVarID, currVarPtr, true);
			WtType tWt = sNode->t->wt;
			if (tWt != -INF) {
				sNode->auxVar = tWt + tFactor;
			} else {
				sNode->auxVar = -INF;
				thenIsZero = true;
			}
			// cout<<"Calling else-inpcf for curNodeID "<<sNode->cnfVarID<<".. "<<std::flush;
			childCnfVarID = inplaceCofactor(sNode->e, sampleCNFVarIDs, currVarPtr);
			WtType eFactor = computeFactor(sampleCNFVarIDs, sNode->cnfVarID, currVarPtr, false);
			WtType eWt = sNode->e->wt;
			if (thenIsZero){
				//if then is zero then we expect else to not be zero. This will be checked in randombits anyway.
				sNode->wt = eWt + eFactor;
			} else {
				if (eWt != -INF) {
					// sNode->wt = sNode->auxVar + eWt * computeFactor(trueCmprsdLvl, cmprsdLvl_e, false, compressedInvPerm);
					WtType op1 = sNode->auxVar, op2 =  eWt + eFactor;
					WtType opMax = fmaxl(op1, op2);
					assert(op1!=-INF); assert(op2!=-INF); assert(opMax!=-INF);
					sNode->wt = log10l(exp10l(op1 - opMax) + exp10l(op2 - opMax)) + opMax;
					assert(sNode->wt!=-INF);
				} else {
					sNode->wt = sNode->auxVar;
				}
			}
		#endif
		ret_CnfVarID = sNode->cnfVarID;
	}
	else {
		//dont compute factor for non-sample-var-nodes. It will be handled in samplevarnodes call to computefactor
		if (t->get(sNode->cnfVarID)){ //var assigned true;
			ret_CnfVarID = inplaceCofactor(sNode->t, sampleCNFVarIDs, currVarPtr);
			sNode->wt = sNode->t->wt;
		} else{//->var assigned false
			ret_CnfVarID = inplaceCofactor(sNode->e, sampleCNFVarIDs, currVarPtr);
			sNode->wt = sNode->e->wt;
		}
		sNode->auxVar = (long) ret_CnfVarID;
	}
	// cout<<"before "<<sNode->cmprsdLvl;
	sNode->cnfVarID = -sNode->cnfVarID;// mark as visited
	// cout<<"after "<<sNode->cmprsdLvl<<"\n";
	// cout<<"Returning from inplacecofactor..\n";
	return ret_CnfVarID;
}

WtType ADDSampler::sampleFromADD(SamplerNode* sNode, vector<Int>& compressedInvPerm, Set<Int>& notSampled){
	if (sNode->t==NULL){
		//snode should be leaf
		// assert(sNode->wt == Cudd_V(sNode->dnode));
		assert(sNode->e==NULL);
		return sNode->wt;
	}
	assert(sNode->cnfVarID < 0);
	if (!t->isSet(-sNode->cnfVarID)){//->var is not newly assigned but needs to be, i.e. its in projectable vars
		// assert(ddVarToCnfVarMap[sNode->dnode->index]==cnfVarIndex);
		#if SAMPLE_NUM_TYPE == 0 || SAMPLE_NUM_TYPE == 2
			bool bit = rb.generateWeightedRandomBit(sNode->auxVar, sNode->wt); // arguments are poswt and totWt
		#else
			bool bit = rb.generateWeightedRandomBit(exp10l(sNode->auxVar), exp10l(sNode->wt)); // arguments are poswt and totWt
			//bool bit = rb.generateWeightedRandomBit(exp10l(sNode->auxVar-sNode->wt), 1);
		#endif
		t->set(bit,-(sNode->cnfVarID));
		assert(notSampled.erase(-(sNode->cnfVarID)));
		numAssigned++;
		return (bit? sampleFromADD(sNode->t, compressedInvPerm, notSampled): sampleFromADD(sNode->e, compressedInvPerm, notSampled));
	} else{
		// assert(ddVarToCnfVarMap[sNode->dnode->index]==trueCNFVarIndex);
		if (t->get(-(sNode->cnfVarID))){//var previously or newly assigned true, check t->get
			return sampleFromADD(sNode->t, compressedInvPerm, notSampled);
		}
		else{//->var previously assigned false
			return sampleFromADD(sNode->e, compressedInvPerm, notSampled);
		}
	}
}

void ADDSampler::inplaceUnCofactor(SamplerNode* sNode){
	// cout<<"Uncofactoring.. "<<std::flush;
	if (sNode->t==NULL){
		//snode should be leaf
		assert(sNode->e==NULL);
		return; // no need to change anything for leaf
	}
	if (sNode->cnfVarID > 0){ // node has been marked unvisited by uncofactor previously
		return;
	}
	//if(cond is true), var is projectable at this node, therefore not previously assigned
	// if (!t->isSet(-(sNode->cnfVarID))){
	// 	inplaceUnCofactor(sNode->t); // will need sn for leaf to access wt!
	// 	inplaceUnCofactor(sNode->e);
	// } else {
	// 	// assert(ddVarToCnfVarMap[sNode->dnode->index]==trueCNFVarIndex);
	// 	if (t->get(-(sNode->cnfVarID))){ //var assigned true;
	// 		inplaceUnCofactor(sNode->t);
	// 	} else{//->var assigned false
	// 		inplaceUnCofactor(sNode->e);
	// 	} 
	// }
	inplaceUnCofactor(sNode->t); // will need sn for leaf to access wt!
	inplaceUnCofactor(sNode->e);
	// cout<<"Before: "<<sNode->cmprsdLvl<<std::flush;
	sNode->cnfVarID = -(sNode->cnfVarID); //mark as unvisited
	// cout<<"After: "<<sNode->cmprsdLvl<<"\n";
}

void ADDSampler::sampleMissing(Set<Int>& notSampled){
	for (auto& cnfvarind: notSampled){
		#if SAMPLE_NUM_TYPE == 0 || SAMPLE_NUM_TYPE == 2
			bool bit = rb.generateWeightedRandomBit(litWts.at(3*cnfvarind), litWts.at(3*cnfvarind+1));
		#elif SAMPLE_NUM_TYPE == 1
			bool bit = rb.generateWeightedRandomBit(exp10l(litWts.at(3*cnfvarind)), exp10l(litWts.at(3*cnfvarind+1)));
		#endif
		t->set(bit,cnfvarind);
		numAssigned++;
	}
}

Asmt& ADDSampler::drawSample(){
	t->clear();
	numAssigned = 0;
	sampleMissing(freeVars);
	//cout<<"Set all freeVars!\n";	
	drawSample_rec(jtRoot);
	assert(numAssigned==nTotalVars);
	return *t;
}

void ADDSampler::drawSample_rec(const JoinNode* jNode){
	if (jNode->projectionVars.size()==0){
		return;
	}
	Dd* a = (Dd*)jNode->dd;
	// cout<<"starting sampling from add\n";
	Set<Int> notSampledVars(jNode->projectionVars);
	a->curVarPtr = a->sampleVarCNFIDs+1;
	// for(Int i = 0; i<a->sampleStructDepth-1; i++){
	// 	cout << a->sampleVarCNFIDs[i+1]<<" ";
	// }
	// cout<<"\n";
	// cout<<"! "<<std::flush;
	inplaceCofactor(rootSNMap.at(jNode), a->sampleVarCNFIDs, &(a->curVarPtr));
	// cout<<"!! "<<std::flush;
	// cout<<"Sampling from "<<jNode->getNodeIndex()+1<<"\n";
	WtType leafVal = sampleFromADD(rootSNMap.at(jNode), a->compressedInvPerm, notSampledVars);
	// cout<<"!!! "<<std::flush;
	sampleMissing(notSampledVars);
	// cout<<"!!!! "<<std::flush;
	inplaceUnCofactor(rootSNMap.at(jNode));
	// cout<<"!!!!! "<<std::flush;
	// Float checkVal = getAsmtVal(a);
	#if SAMPLE_NUM_TYPE != 1
		if (leafVal == 0){
	#else 
		if (leafVal == -INF){
	#endif
		cout<<"While checking sample for add reached leaf with value "<<leafVal<<"\n";
		exit(1);
	} else {
		//cout<<"SampleADD passed!\n";
	}
	for(auto child: jNode->children){
		drawSample_rec(child);
	} 
}

#define TRAVERSE(_cn_,_cl_,_cond_,_defErrMsg_) \
	{while (_cond_){ \
		/*DBG printf("%p ",_cn_);*/ \
		bool _bit_ = t->get(ddVarToCnfVarMap.at(_cl_)); \
		/*DBG printf("%d ",_bit_);*/ \
		_cn_ = _bit_? Cudd_T(_cn_) : Cudd_E(_cn_); \
		_cl_ = _cn_->index; \
	} /*DBGprintf("\n");*/}

double ADDSampler::getAsmtVal(Dd* a){
	DdNode* currNode = a->cuadd.getNode();
	Int currLevel = currNode->index;
	TRAVERSE(currNode,currLevel,!Cudd_IsConstant(currNode),"Error during checkSample traverse");
	//cout<<Cudd_V(currNode)<<" "<<std::flush;
	return (Cudd_V(currNode));
	//return 1;
}

void ADDSampler::writeAsmtToFile(FILE* ofp){
	for(Int i = 1; i<nTotalVars+1; i++){
		fprintf(ofp,"%d ",t->get(i));
	}
	fprintf(ofp,"\n");
	static Int cnt = 0;
	cnt ++;
}

void ADDSampler::sampleAndWriteToFile(FILE* ofp, Int nSamples, bool printProgress){
	Int nTSteps = 10, nT = 1;
	for (Int i =1 ; i<=nSamples; i++){
		drawSample();
		if (printProgress) if (((i-1)%(std::max((nSamples/10),1LL)))==1) util::printComment("Writing trace "+to_string(i)+" out of "+to_string(nSamples)+" to file..");
		writeAsmtToFile(ofp);
		if(i==nT){
			nT *= nTSteps;
		}
	}
}