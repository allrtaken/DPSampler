#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <unordered_map>
#include <utility>

#include "../libraries/cudd/cplusplus/cuddObj.hh"
#include "../libraries/cudd/cudd/cuddInt.h"

#include "../libraries/sylvan/src/sylvan_gmp.h"
#include "../libraries/sylvan/src/sylvan_obj.hpp"

#include "../libraries/cxxopts/include/cxxopts.hpp"

#include "../logic.hh"

using std::unordered_map;
using std::pair;
using std::vector;

namespace SamplerUtils{
//this ADDWrapper is different from addwrapper in aig2dd.cpp. This one is the same as the one in sample_trace
//but renamed to prevent clash with ADD class in cudd. Also some member variable changes. Eventually use namespacce maybe
//changed from struct so making all public. make members private later.
// class ADDWrapper{
// 	public:

// 		ADDWrapper(ADD add_, Int numVars);
// 		ADDWrapper();
// 		ADD add;	
// 		unordered_map<Int, Int> compressedPerm; // stores mapping from ddvarids to compressed-ddvarids/cmprsdlevels
// 		vector<Int> compressedInvPerm; // stores mapping from compressedddvarids/cmprsdlvls to cnfvarid,ddvarid pair
// 		//not explicitly storing mapping from compressed to uncompressed ddvarids. go through cnfvarids for that
// 		Int sampleStructDepth;
// 		bool precompileSampleDAG;
// };

void printPerm(DdManager*);
void printDoubleLine();
// void printComment(const string &message, Int preceedingNewLines, Int followingNewLines, bool commented);
string printTimeTaken(string forWhat, Float timeTaken, Int numNewLines);

// int Cudd_StringSummary(DdManager * dd, DdNode * f, int n, /**< number of variables for minterm computation */
//   string* out /**pointer to string that will contain output*/);
} //end namespace
#endif