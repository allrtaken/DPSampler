#pragma once

/* inclusions =============================================================== */

#include "../libraries/cudd/cplusplus/cuddObj.hh"
#include "../libraries/cudd/cudd/cuddInt.h"

#include "../libraries/sylvan/src/sylvan_gmp.h"
#include "../libraries/sylvan/src/sylvan_obj.hpp"

#include "../libraries/cxxopts/include/cxxopts.hpp"

#include "logic.hh"

/* uses ===================================================================== */

using sylvan::gmp_op_max_CALL;
using sylvan::gmp_op_plus_CALL;
using sylvan::gmp_op_times_CALL;
using sylvan::mtbdd_apply_CALL;
using sylvan::mtbdd_fprintdot_nc;
using sylvan::mtbdd_getdouble;
using sylvan::mtbdd_getvalue;
using sylvan::mtbdd_gmp;
using sylvan::mtbdd_makenode;
using sylvan::Mtbdd;

using cxxopts::value;

/* consts =================================================================== */

const Float MEGA = 1e6; // same as countAntom (1 MB = 1e6 B)

class Dd { // wrapper for CUDD and Sylvan
public:
  ADD cuadd; // CUDD
  Mtbdd mtbdd; // Sylvan

  /*Sampler*/
  std::unordered_map<Int, Int> compressedPerm; // stores mapping from ddvarids to compressed-ddvarids/cmprsdlevels
  vector<Int> compressedInvPerm; // stores mapping from compressedddvarids/cmprsdlvls to cnfvarid,ddvarid pair
  //not explicitly storing mapping from compressed to uncompressed ddvarids. go through cnfvarids for that
  Int* sampleVarCNFIDs;
  Int* curVarPtr;
  Int sampleStructDepth;
  bool precompileSampleDAG;


  Dd(const ADD& cuadd); // CUDD
  Dd(const Mtbdd& mtbdd); // SYLVAN
  Dd(const Dd& dd);

  static const Cudd* newMgr(Float mem, Int threadIndex); // CUDD
  static Dd getConstDd(const Number& n, const Cudd* mgr);
  static Dd getZeroDd(const Cudd* mgr);
  static Dd getOneDd(const Cudd* mgr);
  static Dd getVarDd(Int ddVar, bool val, const Cudd* mgr);
  size_t countNodes() const;
  bool operator<(const Dd& rightDd) const; // *this < rightDd (top of priotity queue is rightmost element)
  Number extractConst() const;
  Dd getComposition(Int ddVar, bool val, const Cudd* mgr) const; // restricts *this to ddVar=val
  Dd getProduct(const Dd& dd) const;
  Dd getSum(const Dd& dd) const;
  Dd getMax(const Dd& dd) const; // real max (not 0-1 max)
  Set<Int> getSupport() const;
  Dd getAbstraction(
    Int ddVar,
    const vector<Int>& ddVarToCnfVarMap,
    const Map<Int, Number>& literalWeights,
    const Assignment& assignment,
    bool additive, // ? getSum : getMax
    const Cudd* mgr
  ) const;
  void writeDotFile(const Cudd* mgr, string dotFileDir = "./") const;
  static void writeInfoFile(const Cudd* mgr, string filePath);
};
