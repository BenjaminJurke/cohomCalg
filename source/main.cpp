////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  main.cpp                                                          +-----------------------+   //
//  ========                                                          | APPLICATION MAIN FILE |   //
//                                                                    +-----------------------+   //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 29.03.2010  File created as main.cpp                                                  //
//                      Contains a fast implementation of the line bundle cohomology algorithm    //
//                      presented in the paper arXiv:1003.5217 [hep-th]                           //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include "main.h"
#include "secondarycohom.h"
#include "rationals.h"
#include "iohandler.h"
#include "platform.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////


#define APP_NAME      "cohomCalg"
#define APP_VERSION   "0.31c"
#define APP_AUTHOR    "author: Benjamin Jurke (mail@benjaminjurke.net)" << std::endl << "    Based on the algorithm presented in arXiv:1003.5217"
#define APP_PLATFORM  APP_TARGET_OS << " " << APP_ARCH


////////////////////////////////////////////////////////////////////////////////////////////////////

vector<char> CBits::cBitsIn16Bits;

void CBits::InitBitCounter()
{
    /* This function initializes the bit counter buffer, by counting the number of bits in all 16-bit
       variables. (Gosh, who remembers the time, when this was ALL of your available memory...?) */

    // Resize the bit buffer
    cBitsIn16Bits.resize((0x1u << 16), 0);

    // And count the bits for all 16-bit variables
    for (unsigned int i=0; i<(0x1u<<16); i++)
    {
       unsigned int count = 0;
       unsigned int n=i;
       while (n)
       {
          count += n & 0x1u;
          n >>= 1;
       }
       cBitsIn16Bits[i] = (char) count;
    }
}

string CBits::IntToBinary(uint64_t x, size_t digits)
{
    /* This functions produces a string containing the lower-end bits of a 64-bit variable up to
       'digits', which are packed into 8-bit packs to ease readability. */

    if ((digits > 0) && (digits <= 64))
    {
        // If we are in range, compute the total length of the output string
        int32_t offset = digits % 8;
        size_t stringlen = digits+digits/8 + (offset == 0 ? -1 : 0);
        string tmp(stringlen, '0');
        uint64_t mask = 0x1ull << (digits-1);

        // Now slide the bit mask over the bits and write a 1 where necessary
        for (int32_t cnt=1,stringpos=0; cnt<=(int32_t) digits; cnt++, stringpos++)
        {
            if (x & mask)
                tmp[stringpos] = '1';

            x <<= 1;
            if(((cnt-offset) % 8 == 0) && (cnt < (int32_t) digits))
                tmp[++stringpos] = ' ';
        }
        return tmp;
    }
    else
    {
        return "";
    }
}

string CBits::IntToBinary(uint64_t x)
{
    /* Returns the bit string containing either 64, 32 or 16 bits depending on numerical value. */

    if (x >= (0x1ull << 32))
        return IntToBinary(x, 64);
    if (x >= (0x1ull << 16))
        return IntToBinary(x, 32);
    return IntToBinary(x, 16);
}


string SecondsToTime(clock_t secs)
{
    /* Converts an integral number of seconds properly into days/hours/minutes/seconds. */

    clock_t days = secs / (60*60*24);
    secs %= 60*60*24;

    clock_t hours = secs / (60*60);
    secs %= 60*60;

    clock_t mins = secs / 60;
    secs %= 60;

    string out;
    char buf[64];
    if (days != 0)  { if (days==1)  out += "1 day ";  else { safe_sprintf(buf, sizeof(buf), "%d days ", (int) days); out += buf; } }
    if (hours != 0) { if (hours==1) out += "1 hour "; else { safe_sprintf(buf, sizeof(buf), "%d hours ", (int) hours); out += buf; } }
    if (mins != 0)  { if (mins==1)  out += "1 min ";  else { safe_sprintf(buf, sizeof(buf), "%d mins ", (int) mins); out += buf; } }
    if ((secs != 0) || out.empty())
                    { if (secs==1)  out += "1 sec ";  else { safe_sprintf(buf, sizeof(buf), "%d secs ", (int) secs); out += buf; } }

    out.erase(out.end()-1);
    return out;
}

string SecondsToTime(double secs)
{
    /* Converts a floating-point number of seconds properly into days/hours/minutes/seconds while
       keeping two digits during the first minute. */

    char buf[256];
    if (secs <= 60.9)
    {
        if (secs == 1.0)
            return "1 second";
        else
            safe_sprintf(buf, sizeof(buf), "%.2f seconds", secs);
    }
    else
        safe_sprintf(buf, sizeof(buf), "%.2f seconds = %s", secs, SecondsToTime((clock_t) secs).c_str());

    return buf;
}

string BytesToReadableSize(uint64_t bytes)
{
    /* Converts a number of bytes to the proper extension string. */

    char buf[32];
    if (bytes >= (1024ull*1024ull*1024ull*1024ull*1024ull))
        safe_sprintf(buf, sizeof(buf), "%.2f PiB", (double) bytes / (1024ull*1024ull*1024ull*1024ull*1024ull));
    else if (bytes >= (1024ull*1024ull*1024ull*1024ull))
        safe_sprintf(buf, sizeof(buf), "%.2f TiB", (double) bytes / (1024ull*1024ull*1024ull*1024ull));
    else if (bytes >= (1024ull*1024ull*1024ull))
        safe_sprintf(buf, sizeof(buf), "%.2f GiB", (double) bytes / (1024ull*1024ull*1024ull));
    else if (bytes >= (1024*1024))
        safe_sprintf(buf, sizeof(buf), "%.2f MiB", (double) bytes / (1024*1024));
    else if (bytes >= 10000)
        safe_sprintf(buf, sizeof(buf), "%.2f KiB", (double) bytes / (1024));
    else if (bytes > 1)
        safe_sprintf(buf, sizeof(buf), "%d bytes", (int) bytes);
    else
        return "1 byte";

    return buf;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


string CCmdLineArguments::strCommand;
string CCmdLineArguments::strInputFileName;
int    CCmdLineArguments::iVerboseLevel = 0;
bool   CCmdLineArguments::bShowTiming = false;
bool   CCmdLineArguments::bShowBits = false;
bool   CCmdLineArguments::bCheckSerre = false;
bool   CCmdLineArguments::bMonomReduction = true;
bool   CCmdLineArguments::bUseMonomFile = true;
bool   CCmdLineArguments::bIntegratedRun = false;
string CCmdLineArguments::strAppendInput;
bool   CCmdLineArguments::bMathematicaOutput = false;
string CCmdLineArguments::strMonomialFileName;
bool   CCmdLineArguments::bMaxVertices = false;
int    CCmdLineArguments::iMaxVertices = -1;
bool   CCmdLineArguments::bMaxSRgens = false;
int    CCmdLineArguments::iMaxSRgens = -1;
bool   CCmdLineArguments::bMaxCohoms = false;
int    CCmdLineArguments::iMaxCohoms = -1;

void CCmdLineArguments::Clear()
{
    strCommand.clear();
    strInputFileName.clear();
    iVerboseLevel = 0;
    bShowTiming = false;
    bCheckSerre = false;
    bMonomReduction = true;
    bUseMonomFile = true;
    strAppendInput.clear();
    bMathematicaOutput = false;
    strMonomialFileName.clear();
}

void CCmdLineArguments::PrintHelp()
{
    /* This function shows the help output. */

    CONSOLE_MSG_OUT("Syntax:  cohomcalg [--option1] [--option2] ... InputFileName [> OutputFileName]");
    CONSOLE_MSG_OUT("");
    CONSOLE_MSG_OUT("Command line options:");
    CONSOLE_MSG_OUT("  --in=\"...\"      Treats the text between parentheses like additional");
    CONSOLE_MSG_OUT("                  input data, i.e. like appended content of the input file.");
    CONSOLE_MSG_OUT("  --nomonomfile   Prohibits usage and generation of monomial file.");
    CONSOLE_MSG_OUT("  --checkserre    Computes the Serre dual cohomology for comparison.");
    CONSOLE_MSG_OUT("  --noreduction   Deactivates the Serre self-duality reduction for ambiguous");
    CONSOLE_MSG_OUT("                  monomials (may increase computation time dramatically!)");
    CONSOLE_MSG_OUT("  --hideinput     Does not print the input data read from the input file.");
    CONSOLE_MSG_OUT("  --showtime      Shows timing stats even for application runs < 1 second.");
    CONSOLE_MSG_OUT("  --showbits      Prints the internally used bitmasks for debug output.");
    CONSOLE_MSG_OUT("  --mathematica   Gives formattet output for the Mathematica script version.");
	CONSOLE_MSG_OUT("  --integrated    Produces minimalistic output for application integration.");
    CONSOLE_MSG_OUT("  --verboseN      Provides debug output at level N=1..6, e.g. --verbose3.");
	CONSOLE_MSG_OUT("  --maxX=N        Defines limits for X=verts,srgens,cohoms to N.");
    CONSOLE_MSG_OUT("");
    CONSOLE_MSG_OUT("The order of options is irrelevant, but later commands always overwrite");
    CONSOLE_MSG_OUT("prior options, e.g. from '--verbose4 --verbose1' only the later one will");
    CONSOLE_MSG_OUT("have any effect.");
    CONSOLE_MSG_OUT("");
    CONSOLE_MSG_OUT("The program automatically tries to add the file extension '.in' to the");
    CONSOLE_MSG_OUT("InputFileName. See the package manual for further information.");
}

bool CCmdLineArguments::ParseCmdLineArguments(int argc, char *argv[])
{
    /* This function parses the C/C++ passed argument lines, which come as an array of strings. Those
       are translated into the static class variables of CCmdLineArguments. */

    Clear();

    // There should at least be one argument
    if (argc < 1)
        return false;

    // The first command is always the executable command
    strCommand = argv[0];

    if (argc < 2)
    {
        // With no option or anything specified, show the help
        PrintHelp();
        return false;
    }

    // Check if someone actually requested the help
    string strCurArg = argv[1];
    if ((strCurArg.compare("--help") == 0) ||
        (strCurArg.compare("-help") == 0) ||
        (strCurArg.compare("/help") == 0) ||
        (strCurArg.compare("--?") == 0) ||
        (strCurArg.compare("-?") == 0) ||
        (strCurArg.compare("/?") == 0))
    {
        PrintHelp();
        return false;
    }

    // The last argument is expected to be a filename, other than than order is irrelevant
    bool bCheckForFilename = true;
    bool bFileRequired = true;
    for (int i=1; i<argc; i++)
    {
        strCurArg = argv[i];

        if (strCurArg.compare(0, 2, "--") == 0)
        {
			/* The verbose levels refer to the level of console output, with a higher number
			   representing more output. Anything below -5 basically means "shut up", i.e. 
			   virtually all output is turned off. */
			if (strCurArg.compare("--integrated") == 0)  { bIntegratedRun = true; continue; }
            if (strCurArg.compare("--hideinput") == 0)   { iVerboseLevel = -1; continue; }
            if (strCurArg.compare("--verbose1") == 0)    { iVerboseLevel =  1; continue; }
            if (strCurArg.compare("--verbose2") == 0)    { iVerboseLevel =  2; continue; }
            if (strCurArg.compare("--verbose3") == 0)    { iVerboseLevel =  3; continue; }
            if (strCurArg.compare("--verbose4") == 0)    { iVerboseLevel =  4; continue; }
            if (strCurArg.compare("--verbose5") == 0)    { iVerboseLevel =  5; continue; }
			if (strCurArg.compare("--verbose6") == 0)    { iVerboseLevel =  6; continue; }
            if (strCurArg.compare("--showbits") == 0)    { bShowBits = true; continue; }
            if (strCurArg.compare("--showtime") == 0)    { bShowTiming = true; continue; }
            if (strCurArg.compare("--checkserre") == 0)  { bCheckSerre = true; continue; }
            if (strCurArg.compare("--noreduction") == 0) { bMonomReduction = false; continue; }
            if (strCurArg.compare("--nomonomfile") == 0) { bUseMonomFile = false; continue; }
            if (strCurArg.compare("--mathematica") == 0) { bMathematicaOutput = true; continue; }
			if (strCurArg.compare(0, 11, "--maxverts=") == 0) { // Read 
				string tmp = argv[i];
				tmp.erase(0, 11);
				int maxverts = atoi(tmp.c_str());
				if ((maxverts < 2) || (maxverts > HARD_MAX_VERTICES)) { ERR_OUT("Invalid value '" << maxverts << "' in --maxverts option, must be 2.." << HARD_MAX_VERTICES); return false; }
				bMaxVertices = true; iMaxVertices = maxverts; 
				continue; 
			}
			if (strCurArg.compare(0, 12, "--maxsrgens=") == 0) { // Read 
				string tmp = argv[i];
				tmp.erase(0, 12);
				int maxsrgens = atoi(tmp.c_str());
				if ((maxsrgens < 1) || (maxsrgens > HARD_MAX_VERTICES)) { ERR_OUT("Invalid value '" << maxsrgens << "' in --maxsrgens option, must be 1.." << HARD_MAX_SRGENS); return false; }
				bMaxSRgens = true; iMaxSRgens = maxsrgens; 
				continue; 
			}
			if (strCurArg.compare(0, 12, "--maxcohoms=") == 0) { // Read 
				string tmp = argv[i];
				tmp.erase(0, 12);
				int maxcohoms = atoi(tmp.c_str());
				if ((maxcohoms < 1) || (maxcohoms > HARD_MAX_VERTICES)) { ERR_OUT("Invalid value '" << maxcohoms << "' in --maxcohoms option, must be 1.." << HARD_MAX_COHOMS); return false; }
				bMaxCohoms = true; iMaxCohoms = maxcohoms; 
				continue; 
			}
            if (strCurArg.compare(0, 5, "--in=") == 0) { // Read in input to append
                strAppendInput = argv[i];
                strAppendInput.erase(0, 5);
                bFileRequired = false;
                continue;
            }
        }
        else
        {
            if ((i == (argc-1)) && bCheckForFilename)
            {
                // The last argument is supposed to be the input file
                string filename = strCurArg;

                // Test if the input file actually exits
                ifstream iFileCheck(filename.c_str());
                if (!iFileCheck)
                {
                    filename += ".in";
                    iFileCheck.open(filename.c_str());
                    if (!iFileCheck)
                    {
                        if (bFileRequired)
                        {
                            ERR_OUT("Could not open file '" << strCurArg << "' or '" << filename << "'.");
                            return false;
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
                iFileCheck.close();

                // Store the filename and the (default) intermediate monomial file name
                strInputFileName = filename;
                strMonomialFileName = strInputFileName + ".monoms";
                continue;
            }
        }

        ERR_OUT("Command line argument " << i << " ('" << strCurArg << "') could not be recognized.");
        CONSOLE_MSG_OUT("Use option '--help' to show a list of available command line options.");
        return false;
    }

	if (bIntegratedRun)
		iVerboseLevel = -9;

    if (bFileRequired && strInputFileName.empty())
    {
        ERR_OUT("No input file specified.");
        return false;
    }

    return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



typedef struct
{
    clock_t t01ApplicationStart;
    clock_t t10SecondarySeqsStart;
    clock_t t20SeqExactnessStart;
    clock_t t30RationalsCountStart;
    clock_t t90ComputationsDone;
    clock_t t99AllDone;
    double  dResolution;
} TimingStats;

////////////////////////////////////////
    int main(int argc, char *argv[])    //            ###  APPLICATION MAIN ROUTINE  ###
////////////////////////////////////////
{
    /* Well, the main function probably does not require an explanation, but ok: Basically, the
       main routine tells the compiler where the code entry point of the resulting binary is to
       be found and where to start executing the commands. This also defines the lowest level of
       the stack for the running application. So this is, where it all begins.... */

    // Start with aquiring some timing information
    TimingStats ts;
    ts.t01ApplicationStart = clock();
    ts.dResolution = 1.0 / CLOCKS_PER_SEC;

    // Output the application information
    stringstream ss;
    ss << endl
       << "    " << APP_NAME << " v" << APP_VERSION << endl
       << "    (compiled for " << APP_PLATFORM << ")" << endl
       << "    " << APP_AUTHOR << endl
       << endl;
    string strAppInfo = ss.str();

    CONSOLE_MSG_OUT(strAppInfo);

    // STEP 1: Analyze the command line arguments
    if (!CCmdLineArguments::ParseCmdLineArguments(argc, argv))
    {
        CONSOLE_MSG_OUT("");
		if (CCmdLineArguments::GetVerboseLevel() < -5)
			MSG_OUT_NOENDL("{False,\"Invalid command line parameters\"}");
        return 1;
    }

	if (CCmdLineArguments::GetVerboseLevel() >= -5)
		MSG_OUT(strAppInfo);


    // STEP 2: Do all pre-computation initializations
    CBits::InitBitCounter();


    //////////////////////////


    // STEP 3: Read in the input file and print the input data
    CInternalData id;
    if (!id.ReadAndParseInputFile(CCmdLineArguments::GetInputFileName(), CCmdLineArguments::GetAppendInput()))
	{
		if (CCmdLineArguments::GetVerboseLevel() < -5)
			MSG_OUT_NOENDL("{False,\"Invalid input data\"}");
        return 1;
	}

    if (CCmdLineArguments::GetVerboseLevel() >= 0)
    {
        MSG_OUT("Input data:");
        MSG_OUT("===========");
        id.PrintInternalData();
    }

	if (CCmdLineArguments::GetVerboseLevel() >= -5)
	{
		MSG_OUT("");
		MSG_OUT("");
		MSG_OUT("");
	}

    //////////////////////////


    // STEP 4: Initialize the secondary sequences engine
    CSecondaryCohomology seccohom;
    if (!seccohom.Init2ndSequences(id))
	{
		if (CCmdLineArguments::GetVerboseLevel() < -5)
			MSG_OUT_NOENDL("{False,\"Could not compute secondary sequences\"}");
        return 1;
	}

    ts.t10SecondarySeqsStart = clock();

    // Check if we should use a potentially available monomial file for speedup, if so skip steps 5-6
    string strMonomialFileName = CCmdLineArguments::GetMonomialFileName();
    if (strMonomialFileName.empty())
        CCmdLineArguments::SetUseMonomFile(false);

    bool bComputeMonoms = true;
    if (CCmdLineArguments::GetUseMonomFile())
    {
        if (seccohom.GetMonomialsList().ReadMonomialsFile(id, strMonomialFileName))
            bComputeMonoms = false;
        else
        {
            CONSOLE_MSG_OUT("Could not open/read in the monomial file '" << strMonomialFileName << "'.");
            CONSOLE_MSG_OUT("");
        }
    }
    else
    {
        CONSOLE_MSG_OUT("Usage and generation of intermediate monomial files deactivated.");
        CONSOLE_MSG_OUT("");
    }

    if (bComputeMonoms)
    {
        // STEP 5: Start the computation of the secondary sequences
        CONSOLE_MSG_OUT("Starting computation of secondary sequences...");
		if (CCmdLineArguments::GetVerboseLevel() >= 6)
			MSG_OUT("   - Generating secondary sequences...");

        if (!seccohom.TraverseSRpowerset(id))
		{
			if (CCmdLineArguments::GetVerboseLevel() < -5)
				MSG_OUT_NOENDL("{False,\"Could not traverse SR ideal generator powerset\"}");
            return 1;
		}

        if (CCmdLineArguments::GetVerboseLevel() >= 4)
        {
            MSG_OUT("Verbose Level 4: Full list of monomials with secondary sequences:");
            MSG_OUT("-----------------------------------------------------------------");
            seccohom.PrintMonomMap(id, false);
            MSG_OUT("");
            MSG_OUT("");
        }

        // STEP 6a: Compute the secondary cohomology / monomial factors from the sequences
		if (CCmdLineArguments::GetVerboseLevel() >= 6)
			MSG_OUT("   - Eliminating trivial sequences...");

        ts.t20SeqExactnessStart = clock();
        if (!seccohom.Compute2ndCohomFromTrivialSequences())
		{
			if (CCmdLineArguments::GetVerboseLevel() < -5)
				MSG_OUT_NOENDL("{False,\"Could not compute secondary cohomology of trivial sequences\"}");
            return 1;
		}

        if (CCmdLineArguments::GetVerboseLevel() >= 3)
        {
            MSG_OUT("Verbose Level 3: Reduced list of monomials with non-trivial secondary sequences:");
            MSG_OUT("--------------------------------------------------------------------------------");
            seccohom.PrintMonomMap(id, false);
            MSG_OUT("");
            MSG_OUT("");
        }

        // STEP 6b: Compute the remaining secondary/remnant cohomology elements
		if (CCmdLineArguments::GetVerboseLevel() >= 6)
			MSG_OUT("   - Computing cohomology from sequences...");
        if (!seccohom.Compute2ndCohomFromSequences(id))
		{
			if (CCmdLineArguments::GetVerboseLevel() < -5)
				MSG_OUT_NOENDL("{False,\"Could not compute secondary cohomology\"}");
            return 1;
		}

        // If we are using intermediate monomial files, it is time for storing the computed data to file
        if (CCmdLineArguments::GetUseMonomFile())
        {
            if (!seccohom.GetMonomialsList().WriteMonomialsFile(id, strMonomialFileName))
            {
                ERR_OUT("Could not write the computed monomials to the file '" << strMonomialFileName << "'. Continuing...");
            }
        }
    }
    else
    {
        CONSOLE_MSG_OUT("Using monomial file '" << strMonomialFileName << "'. Delete or rename this file");
        CONSOLE_MSG_OUT("if you don't want to use it, or use '--nomonomfile' option.");
        CONSOLE_MSG_OUT("WARNING: Do not change the vertices (order) or Stanley-Reisner ideal in the");
        CONSOLE_MSG_OUT("input file without deleting the monomial file - otherwise crashes!");
        CONSOLE_MSG_OUT("");
    }

    // STEP 7: Resolve the ambiguous monomials via the "unique complement dualization" technique (unless turned off)
    if (CCmdLineArguments::GetMonomReduction())
    {
        if (CCmdLineArguments::GetVerboseLevel() >= 2)
        {
            MSG_OUT("Verbose Level 2: Preliminary list of contributing monomials with factors (pre-Serre-reduction):");
            MSG_OUT("-----------------------------------------------------------------------------------------------");
            seccohom.GetMonomialsList().PrintMonomialList(id, true, false, false);
            MSG_OUT("");
            MSG_OUT("");
        }

        if (!seccohom.Perform2ndCohomSerreReduction(id))
		{
			if (CCmdLineArguments::GetVerboseLevel() < -5)
				MSG_OUT_NOENDL("{False,\"Could not perform Serre duality reduction of contributions\"}");
            return 1;
		}
		if (CCmdLineArguments::GetVerboseLevel() >= -5)
		{
			MSG_OUT("");
			MSG_OUT("");
		}
    }
    else
    {
        CONSOLE_MSG_OUT("Serre self-duality reduction for ambiguous monomials deactivated.");
    }

    if (CCmdLineArguments::GetVerboseLevel() >= 2)
    {
        MSG_OUT("Verbose Level 2: Final list of contributing monomials with factors:");
        MSG_OUT("-------------------------------------------------------------------");
        seccohom.GetMonomialsList().PrintMonomialList(id, true, false, false);
        MSG_OUT("");
        MSG_OUT("");
    }

    CONSOLE_MSG_OUT("Computation of secondary cohomologies and contributions complete.");

    //////////////////////////

    // STEP 8: Count the number of rational functions for alle monomials
    ts.t30RationalsCountStart = clock();

    vector<CCohomology> cohom;
    if (!CRationals::ComputeCohomologies(id, seccohom.GetMonomialsList(), cohom))
    {
		if (CCmdLineArguments::GetVerboseLevel() < -5)
			MSG_OUT_NOENDL("{False,\"Internal error while computing the cohomology\"}");
		else
			MSG_OUT("Internal error while computing the cohomology");
        return 1;
    }

    ts.t90ComputationsDone = clock();


    if (CCmdLineArguments::GetVerboseLevel() >= 1)
    {
        MSG_OUT("Verbose Level 1: Final list of contributing monomials with factors and rationals:");
        MSG_OUT("---------------------------------------------------------------------------------");
        for (unsigned int i=0; i<cohom.size(); i++)
        {
            cohom[i].PrintFullCohomologyMonomialMap(id, true);
            MSG_OUT("");
        }
        MSG_OUT("");
    }

    if (CCmdLineArguments::GetVerboseLevel() >= -5)
	    MSG_OUT("");

    //////////////////////////

    // STEP 9: Provide the final output

    if (CCmdLineArguments::GetMathematicaOutput())
    {
        MSG_OUT("Preformatted output for Mathematica script version of the algorithm implementation:");
        MSG_OUT("-----------------------------------------------------------------------------------");

        string tmp;
        MSG_OUT(id.PrintInternalDataAsMathematicaScriptInput());
        CCohomology::GetMathematicaCohomologiesList(cohom, tmp);
        MSG_OUT(tmp);
        MSG_OUT("");
        MSG_OUT("");
    }

	if (!CCmdLineArguments::GetIntegratedRun())
		CCohomology::PrintCohomologies(cohom);
	else
		CCohomology::SummarizeCohomologies(cohom);

    ts.t99AllDone = clock();


    // STEP 10: Output timing statistics
    double totaltime = (ts.t99AllDone-ts.t01ApplicationStart)*ts.dResolution;
    if (((totaltime > 1.0) || (CCmdLineArguments::GetShowTiming())) && (!CCmdLineArguments::GetIntegratedRun()))
    {
        double secondarycohoms = (ts.t30RationalsCountStart - ts.t10SecondarySeqsStart)*ts.dResolution;
        double rationalscount  = (ts.t90ComputationsDone - ts.t30RationalsCountStart)*ts.dResolution;
        MSG_OUT("");
        MSG_OUT("");
        MSG_OUT("Application run took " << SecondsToTime(totaltime) << ", more precisely");
        MSG_OUT("    " << SecondsToTime(secondarycohoms) << " for the computation of the secondary cohomology");
        MSG_OUT("    " << SecondsToTime(rationalscount)  << " for the counting of rational functions");
    }

    CONSOLE_MSG_OUT("");
    CONSOLE_MSG_OUT("    All done. Program run successfully completed.");
    CONSOLE_MSG_OUT("");

    return 0;
}
