////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  rationals.h                                                                                   //
//  ===========                                                                                   //
//                                                                                                //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 08.04.2010  File created as rationals.h                                               //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef INC_RATIONALS_H
#define INC_RATIONALS_H


#include <string>
#include <map>
#include <vector>
#include <thread>
#include "iohandler.h"
#include "secondarycohom.h"


////////////////////////////////////////////////////////////////////////////////////////////////////

// In this file we take the list of monomials and secondary cohomology factors determined earlier
// and count the number of rational functions with the appropiate variables in the numerator and
// denominator. More precisely, this problem can be reformulated as a underdetermined linear system
// of equations, which is solved over the positive integers. Therefore the counting of rational 
// functions is equivalent to counting the number of integer lattice points inside a corresponding
// polyhedron determined by those constraints. We rely on the freely available PolyLib library for
// the task of counting those integers, which makes use of the theory of Ehrhard polynomials.
// 
// As for the data structures used in this class: The class CCohomology is the "result" class, which
// contains the requested charges of the divisor determining the line bundle O(D) again and of 
// course the computed dimensions of the cohomology classes H^i(X;O(D)).
//
// The actual computations are carried out in the CRationals class, which basically serves as a
// wrapper for the PolyLib and translates our data format appropiately.


class CCohomology
{
  friend class CRationals;

  private:
    // The output data
    std::vector<int32_t>                BundleGLSMch;     // determines divisor D of line bundle O(D)
    std::vector<std::vector<uint32_t> > CohomologyDims;   // contains (candidate) dimensions of H^i(X;O(D))
	std::vector<std::string>            ContributingDenoms; // contains a vector of contributing denominator monomials

    CMonomialsList monoms;                  // contains all monomials and their contributions to the results

  public:
    // Output functions
    static std::string GetCohomologyGroupString(const std::vector<int32_t> &BundleGLSMcharges);
    static std::string GetGroupDimensionsString(const std::vector<uint32_t> &CohomologyDimensions);
	static std::string GetGroupDimensionsStringNoPadding(const std::vector<uint32_t> &CohomologyDimensions);
    static std::string GetCohomologyString(const std::vector<int32_t> &BundleGLSMcharges, const std::vector<uint32_t> &CohomologyDimensions);
    std::string GetCohomologyString() const;

    void PrintFullCohomologyMonomialMap(const CInternalData &id, bool ShortList) const;
    static void PrintCohomologies(const std::vector<CCohomology> &cohomologies);
	static void SummarizeCohomologies(const std::vector<CCohomology> &cohomologies);
    static void GetMathematicaCohomologiesList(const std::vector<CCohomology> &cohomologies, std::string &out);

	// Data retrieval
    inline bool IsAmbiguous() const     { return (CohomologyDims.size() != 1); };
    inline bool IsNotDetermined() const { return (CohomologyDims.size() < 1); };
};


class CRationals
{
  private:
	static std::vector< std::pair<std::thread *, void *> > WorkersList;       // contains the list of worker threads;
	static bool WaitForAllWorkersFinish();

  private:
    static void PrintPolyLibConstraintMatrix(const CInternalData &id, void *m);
    static bool CountRationalFctns(const CInternalData &id, uint64_t monomial, const i32vec64 &TargetDivisor, uint32_t &out);

    static bool ComputeCohomology(const CInternalData &id, CMonomialsList &ml, const i32vec64 &TargetDivisor, std::vector<uint32_t> &out_cohomology, std::vector<std::string> &out_contribcohoms);
    static bool ComputeCohomologyAndSerreDual(const CInternalData &id, CMonomialsList &ml, const i32vec64 &TargetDivisor, std::vector<uint32_t> &out_cohomology, std::vector<std::string> &out_contribcohoms);

  public:
    static bool ComputeCohomologies(const CInternalData &id, const CMonomialsList &ml, std::vector<CCohomology> &out_cohomologies);
};


#endif