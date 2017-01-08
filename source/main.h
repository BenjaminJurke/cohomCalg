////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  main.h                                                            +-----------------------+   //
//  ======                                                            | APPLICATION MAIN FILE |   //
//                                                                    +-----------------------+   //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 31.03.2010  File created as main.h                                                    //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef INC_MAIN_H
#define INC_MAIN_H


#include <stdint.h>
#include <string>
#include <iostream>
#include <vector>
#include <ctime>


////////////////////////////////////////////////////////////////////////////////////////////////////

// This is the main header file, which is included by all other source files. It contains a couple of
// "helper" functions for the conversion of seconds to readable time, for counting and properly printing
// of bitmasks. Furthermore, it also contains the CCmdLineArguments class, which allows access to the
// command line parameters of the program. This header file also declares the primary output functions
// (or rather macros) which are used throughout the program.


// A couple of hard internal limits
#define HARD_MAX_VERTICES 64
#define HARD_MAX_SRGENS 64
#define HARD_MAX_COHOMS 1024


// Yeah, I know, macros are evil... spare me the lecture
#define X_COUNT_BITS16(___n_)   (CBits::cBitsIn16Bits[___n_])
#define X_COUNT_BITS32(___n_)   (CBits::cBitsIn16Bits[___n_ & 0xffffu] + CBits::cBitsIn16Bits[(___n_ >> 16) & 0xffffu])
#define X_COUNT_BITS64(___n_)   (CBits::cBitsIn16Bits[___n_ & 0xffffu] + CBits::cBitsIn16Bits[(___n_ >> 16) & 0xffffu] + CBits::cBitsIn16Bits[(___n_ >> 32) & 0xffffu] + CBits::cBitsIn16Bits[(___n_ >> 48) & 0xffffu])

class CBits
{
  private:
    // The bit buffer
    static std::vector<char> cBitsIn16Bits;

  public:
    // Initialization
    static void InitBitCounter();

    // Conversion and bit counting
    static std::string IntToBinary(uint64_t x, size_t digits);
    static std::string IntToBinary(uint64_t x);
    static inline unsigned int CountBits(uint16_t x)  { return X_COUNT_BITS16(x); }
    static inline unsigned int CountBits(uint32_t x)  { return X_COUNT_BITS32(x); }
    static inline unsigned int CountBits(uint64_t x)  { return X_COUNT_BITS64(x); }
};


// A couple of helper functions
std::string SecondsToTime(clock_t secs);
std::string SecondsToTime(double secs);
std::string BytesToReadableSize(uint64_t bytes);


////////////////////////////////////////////////////////////////////////////////////////////////////


class CCmdLineArguments
{
  private:
    // Command line parameter values
    static std::string strCommand;
    static std::string strInputFileName;
    static int         iVerboseLevel;
    static bool        bShowTiming;
    static bool        bShowBits;
    static bool        bCheckSerre;
    static bool        bMonomReduction;
    static bool        bUseMonomFile;
	static bool        bIntegratedRun;
    static std::string strAppendInput;
    static bool        bMathematicaOutput;
    static std::string strMonomialFileName;
	static bool        bMaxVertices;
	static int         iMaxVertices;
	static bool        bMaxSRgens;
	static int         iMaxSRgens;
	static bool        bMaxCohoms;
	static int         iMaxCohoms;

  private:
    static void Clear();
    static void PrintHelp();

  public:
    // Primary parsing
    static bool ParseCmdLineArguments(int argc, char **argv);

    // Almost all command line variables are read-only
    static inline std::string GetInputFileName()     { return strInputFileName; };
    static inline int         GetVerboseLevel()      { return iVerboseLevel; };
    static inline bool        GetShowTiming()        { return bShowTiming; };
    static inline bool        GetShowBits()          { return bShowBits; };
    static inline bool        GetCheckSerre()        { return bCheckSerre; };
    static inline bool        GetMonomReduction()    { return bMonomReduction; };
    static inline bool        GetUseMonomFile()      { return bUseMonomFile; };
    static inline std::string GetAppendInput()       { return strAppendInput; };
    static inline bool        GetMathematicaOutput() { return bMathematicaOutput; };
	static inline bool        GetIntegratedRun()     { return bIntegratedRun; };
    static inline std::string GetMonomialFileName()  { return strMonomialFileName; };
	static inline int         GetMaxVertices()       { if (bMaxVertices) return iMaxVertices; else return HARD_MAX_VERTICES; }
	static inline int         GetMaxSRgens()         { if (bMaxSRgens) return iMaxSRgens; else return HARD_MAX_SRGENS; }
	static inline int         GetMaxCohoms()         { if (bMaxCohoms) return iMaxCohoms; else return HARD_MAX_COHOMS; }
    
    // Only the monomial filename and usage can be changed
    static inline void        SetMonomialFileName(std::string &filename) { strMonomialFileName = filename; };
    static inline void        SetUseMonomFile(bool status)               { bUseMonomFile = status; };
};


////////////////////////////////////////////////////////////////////////////////////////////////////

// The primary output macros
#define CONSOLE_OUT(msg)      std::cerr << msg
#define CONSOLE_MSG_OUT(msg)  CONSOLE_OUT(msg << std::endl)
#define STATUS_OUT(msg)       CONSOLE_OUT("STATUS: " << msg << "\r")
#define ERR_OUT_PLAIN(msg)    CONSOLE_OUT(msg << std::endl)
#define ERR_OUT(errmsg)       ERR_OUT_PLAIN("ERROR: " << errmsg)
#define WARN_OUT(errmsg)      ERR_OUT_PLAIN("WARNING: " << errmsg)
#define MSG_OUT_NOENDL(msg)   std::cout << msg
#define MSG_OUT(msg)          MSG_OUT_NOENDL(msg << std::endl)


#endif
