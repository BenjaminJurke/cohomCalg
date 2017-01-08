////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  iohandler.h                                                                                   //
//  ===========                                                                                   //
//                                                                                                //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 17.04.2010  File created as iohandler.h                                               //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef INC_IOHANDLER_H
#define INC_IOHANDLER_H


#include <stdint.h>
#include <string>
#include <vector>
#include <set>


////////////////////////////////////////////////////////////////////////////////////////////////////

// This file handles the translation from the input file or string to the internal data format using
// first the tokenizer and then some pretty elementary (read: non-sophisticated) parsing.
//
// There are several internal data format and due to the quite limited amount of input data several
// redundant copies are held by the computation classes. However, the CInternalData class provides
// the "common ground" and also contains the output functions to "human readable" format.
//
// Internally, for most of the program we are working with 64-bit integer variables, because the
// usage of bitmasks and bit-wise operators in general is very fast and rather memory efficient.
// The usage of those internal variables also puts some limitations on the input data: only up to
// 63 coordinates/vertices and Stanley-Reisner ideal generator are supported (the top bit is
// reserved!), butat the current state of technology those cases are way out of the league of modern
// computers.


// For speedup we use some fixed-width arrays of 64 integer variables
template <class T> class vec64
{
  public:
    T x[64];

  public:
    inline void Fill(T val) { for (unsigned int i=0; i<64; i++) x[i]=val; }
    inline void Clear()     { Fill(0); };
};

typedef vec64<int32_t>  i32vec64;
typedef vec64<uint32_t> ui32vec64;

typedef vec64<int64_t>  i64vec64;
typedef vec64<uint64_t> ui64vec64;

// For each coordinate/vertex we hold the following data
typedef struct
{
    uint64_t    liVar;
    std::string strName;
    i32vec64    GLSMch;
} InternalCoordData;


class CInternalData
{
  friend class CInputFile;

  private:
    // The input data
    std::vector<InternalCoordData> vCoords;         // The coordinates/vertices
    std::vector<uint64_t>          vSRgens;         // The Stanley-Reisner ideal generators
    std::vector<i32vec64>          vTargetDivisors; // The requested ambient cohomology bundle charges
    size_t                         numGLSMch;

    // For output formatting
    size_t nMaxMonomWidth;

    // Some internal data
    uint64_t liCompleteUnion;

  private:
    void Clear();

  public:
    CInternalData();

    // Output functions
    std::string Int64ToCoordProduct(uint64_t liProduct, std::string strSep = "*", std::string strZeroVal = "1") const;
    std::string Int64ToCoordProductPadded(uint64_t liProduct, std::string strSep = "*", std::string strZeroVal = "1") const;
    std::string Int64ToMonomial(uint64_t liProduct) const;
    std::string Int64ToMonomialPadded(uint64_t liProduct) const;

    void        PrintInternalData();
    std::string PrintInternalDataAsMathematicaScriptInput();

    // Internal data retrival
    inline size_t                                GetDimension() const         { return vCoords.size() - numGLSMch; };
    inline uint64_t                              GetCompleteUnion() const     { return liCompleteUnion; };
    inline size_t                                GetNumGLSMch() const         { return numGLSMch; };
    inline size_t                                GetNumCoordinates() const    { return vCoords.size(); };
    inline size_t                                GetNumSRgens() const         { return vSRgens.size(); };
    inline size_t                                GetNumTargetDivisors() const { return vTargetDivisors.size(); };
    inline const std::vector<InternalCoordData> &GetInternalCoordData() const { return vCoords; };
    inline const std::vector<uint64_t> &         GetSRgens() const            { return vSRgens; };
    inline const std::vector<i32vec64> &         GetTargetDivisors() const    { return vTargetDivisors; };

    void GetCanonicalDivisor(i32vec64 &candiv_out) const;

    // Data input
    bool ReadAndParseInputFile(std::string strFileName, std::string strAppend);
};


////////////////////////////////////////////////////////////////////////////////////////////////////

// The CMonomialsList class contains the data which is ultimately produced by the processes. Encoded
// in 64-bit variables, again, it contains the list of all contributing monomials and their respective
// factor, which was derived from the secondary/remnant cohomology. Note that the monomials are split
// up in two sets: the monomials which have a contribution to an uniquely determined cohomology group
// and the ambiguously contributing monomials, whose actual contribution has to be sorted out later.

typedef struct
{
    uint32_t nGroup;
    uint32_t nFactor;
} CohomContrib;

typedef struct
{
    CohomContrib Cohom;
    uint32_t nRationals;
    uint32_t nRationalsDual;
} UniqueContribData;

typedef struct
{
    std::vector<CohomContrib> vCohoms;
    uint32_t nRationals;
    uint32_t nRationalsDual;
} AmbiguousContribData;

class CMonomialsList
{
  friend class CSecondaryCohomology;
  friend class CRationals;
  friend class CMonomialFile;

  private:
    // Associative monomial contribution data
    std::map<uint64_t, UniqueContribData>    unique_monoms;
    std::map<uint64_t, AmbiguousContribData> ambiguous_monoms;

  private:
    // Adding monomials
    bool AddUniqueContribution(uint64_t liMonomial, uint32_t nCohomGroup, uint32_t nFactor);
    bool AddAmbiguousContribution(uint64_t liMonomial, const std::vector<CohomContrib> &vContributions);

  public:
    void Clear();
    void ClearRationals();

    // Output functions
    void PrintMonomialList(const CInternalData &id, bool bPrintFactors, bool bPrintRationals, bool bShortList) const;

    // Intermediate monomial file I/O
    bool ReadMonomialsFile(const CInternalData &id, std::string strFileName);
    bool WriteMonomialsFile(const CInternalData &id, std::string strFileName) const;
};

bool GetIntSequenceStartEnd(const std::vector<uint64_t> &seq, size_t &min_out, size_t &max_out);

#endif
