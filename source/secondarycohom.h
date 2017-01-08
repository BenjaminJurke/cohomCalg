////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  secondarycohom.h                                                                              //
//  ================                                                                              //
//                                                                                                //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 04.04.2010  File created as secondarycohom.cpp                                        //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef INC_SECONDARYCOHOM_H
#define INC_SECONDARYCOHOM_H


#include <stdint.h>
#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include "iohandler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////

// This first computational part of the program houses the engine that traverses the powerset of the
// Stanley-Reisner ideal generators. Using 64-bit integer bit-masks, which are specified in the
// CInternalData handler, is basically runs a loop through all possible combinations and successively
// builds up the secondary/remnant sequences. After the loop, the cohomology of those sequences is
// determined by considering the alternating sum - which is a quick exactness test - and attempting
// exactness repairs in the range of C^0 to C^dim of the variety, which ultimately determine the
// secondary cohomology. All those steps are handled by the CSecondaryCohomology class.


typedef struct
{
    uint64_t num;
    int64_t  altersum;
    std::vector<uint64_t> c;

  public:
    inline void Clear(size_t len) { num = 0; altersum = 0; c.clear(); if (len>0) c.resize(len, 0); };
} MonomData;

typedef std::map<uint64_t, MonomData> MonomMap;

class CSecondaryCohomology
{
  private:
    // Secondary sequences data
    MonomData insertion_prototype;
    MonomMap  unique_monoms;

    int32_t min_c_deg;
    int32_t max_c_deg;
    uint32_t len_c_deg;

    CMonomialsList monoms;

  private:
    // Internal computations
    bool IsSequenceExact(const std::vector<uint64_t> &seq, size_t &pos_out) const;

    // Internal output functions
    std::string Print2ndSeq(const MonomData &mondat, int32_t min_out, int32_t max_out);
    std::string Print2ndSeq(const MonomData &mondat);
    std::string PrintMonomialWith2ndSeq(const CInternalData &id, MonomMap::const_iterator monom, int min_out, int max_out);

    void Clear();

  public:
    CSecondaryCohomology();

    // The primary SR powerset traverse functions (EXECUTE IN ORDER!)
    bool Init2ndSequences(const CInternalData &id);              // STEP 1
    bool TraverseSRpowerset(const CInternalData &id);            // STEP 2
    bool Compute2ndCohomFromTrivialSequences();                  // STEP 3
    bool Compute2ndCohomFromSequences(const CInternalData &id);  // STEP 4
    bool Perform2ndCohomSerreReduction(const CInternalData &id); // STEP 5

    // Monomial list access
    inline CMonomialsList &GetMonomialsList() { return monoms; };

    // Output function
    void PrintMonomMap(const CInternalData &id, bool shorten_output);
};


#endif
