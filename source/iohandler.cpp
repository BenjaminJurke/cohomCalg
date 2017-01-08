////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  iohandler.cpp                                                                                 //
//  =============                                                                                 //
//                                                                                                //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 17.04.2010  File created as iohandler.cpp                                             //
//                      Parses and analyzes the input file format and converts to the internal    //
//                      fast 64-bit integer formal                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <map>

#include "iohandler.h"
#include "tokenizer.h"
#include "main.h"
#include "platform.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////


#define MAX_INPUT_LINES   10000         // Maximal number of lines read from an input file
#define MAX_GLSM_RANGE    (1024*1024)   // GLSM charges may range from -(val) to +(val)


////////////////////////////////////////////////////////////////////////////////////////////////////

// Besides the class CInternalData defined in the header, there is the class CInputFile which
// entirely handles the interpretation of the input data using the tokenizer, i.e. this class
// basically serves as the "parser" and also handles the first round of checking of the input data.
// Note that this class is entirely restricted to this input file and should never be handles
// externally.

class CInputFile
{
  private:
    // Read data
    map< string, vector<int32_t> > mCoords;         // Coordinates and GLSM charges
    vector< vector<string> >       vSRideal;        // Stanley-Reisner ideal generators
    vector< vector<int32_t> >      vAmbientCohoms;  // Requested ambient cohomologies
    size_t                         numGLSMch;       // Number of GLSM charges (for syntax check)

  private:
    // The internal parsing functions
    bool ParseVertexCommand(CTokenizer &tok);
    bool ParseSRidealCommand(CTokenizer &tok);
    bool ParseAmbientCohomCommand(CTokenizer &tok);
    bool ParseMonomialFileCommand(CTokenizer &tok);
	bool DiscardCommand(CTokenizer &tok);
    bool ParseInput(const char *pInputData);

    void HandleSyntaxError(size_t iCommandNumber, const char *pCmdStart, ptrdiff_t iErrOffset);

    bool TranslateToInternalData(CInternalData &intdata);

  public:
    CInputFile();
    void Clear();

    // Main data handling functions
    static bool ReadInputFileWithoutComments(string strFileName, string &out_contents);
    bool        ReadAndParseInputFile(string strFileName, string strAppend, CInternalData &out);
};


CInputFile::CInputFile()
{
    Clear();
}

void CInputFile::Clear()
{
    mCoords.clear();
    vSRideal.clear();
    vAmbientCohoms.clear();
    numGLSMch = 0;
}

bool CInputFile::ReadInputFileWithoutComments(string strFileName, string &out_contents)
{
    /* This function simply reads in all lines from the text file and cuts of everything after
       a '%' character, which is reserved for comments. The result is a "comment-free" string
       containing the input file.
       Note that only MAX_INPUT_LINES number of lines (see definition at the top of file) are
       retrieved from the file. */

    ifstream ifs(strFileName.c_str());
    if (!ifs.is_open())
    {
        ERR_OUT("Could not open or read the file '" << strFileName << "'.");
        return false;
    }

    out_contents.clear();

    string line;

    unsigned int nLines = 0;

    // Work through the lines and cut of everything after a '%'-character
    while ((getline(ifs, line)) && (nLines++ < MAX_INPUT_LINES))
    {
        size_t comment = line.find('%');
        out_contents.append("\n");
        if (comment != string::npos)
            out_contents.append(line.substr(0,comment));
        else
            out_contents.append(line);
    }

    // Have we reached the legal maximum of lines - which is considered an error!
    if (nLines >= MAX_INPUT_LINES)
    {
        ERR_OUT("The input file '" << strFileName << "' contains more than " << MAX_INPUT_LINES << " lines.");
        return false;
    }

    return true;
}


bool CInputFile::ReadAndParseInputFile(string strFileName, string strAppend, CInternalData &out)
{
    /* As one may provide input data by file, command line or both, this functions takes care
       of putting those two input sources correctly together. The resulting input is then
       subject to the parsing functions and finally translated into the internal data format
       described in the class CInternalData. */

    string input_data;

    if (!strFileName.empty())
    {
        if (!ReadInputFileWithoutComments(strFileName, input_data))
            return false;
    }

    input_data.append(strAppend);

    if (!ParseInput(input_data.c_str()))
        return false;

    if (!TranslateToInternalData(out))
        return false;

    return true;
}


void CInputFile::HandleSyntaxError(size_t iCommandNumber, const char *pCmdStart, ptrdiff_t iErrOffset)
{
    /* In case the parsing functions stumble upon some unexpected input, this function handles
       a properly formatted output of the command line indicating where the error was raised. */

    if (iErrOffset < 0)
        iErrOffset = 0;

    string error_line;
    ptrdiff_t error_line_offset;
    if (iErrOffset > 45)
    {
        // If the error comes up at a late character in the line, shorten
        error_line.assign(pCmdStart, 16);
        error_line.append("  ...  ");
        error_line.append(pCmdStart + iErrOffset - 15, 30);
        error_line_offset = 16 + 7 + 15;
    }
    else
    {
        error_line.assign(pCmdStart, iErrOffset + 15);
        error_line_offset = iErrOffset;
    }

    // Truncate at newline or carriage return characters that may appear after the error
    size_t oddch = error_line.find_first_of("\n\r", error_line_offset);
    if (oddch != string::npos)
        error_line.erase(oddch);

    // Replace all newline, tab or carriage return characters in the final string by spaces
    oddch = error_line.find_first_of("\t\n\r", error_line_offset);
    while (oddch != string::npos)
    {
        error_line[oddch] = ' ';
        oddch = error_line.find_first_of("\t\n\r", error_line_offset);
    }

    // Print the error
    ERR_OUT("Syntax error in input file, command no. " << iCommandNumber << ":");
    ERR_OUT_PLAIN(error_line);
    ERR_OUT_PLAIN(string(error_line_offset, ' ') << '^');
}


bool CInputFile::ParseVertexCommand(CTokenizer &tok)
{
    /* This function parses a 'vertex' command, for which the maximal prototype syntax line is

       vertex z  = (  0,  0,  0, -2, -3) | PIC: F | GLSM: (-4,  1,  2,  0,  0,  0,  0,  0);
                 ********************************
       vertex z  = (  0,  0,  0, -2, -3) | GLSM: (-4,  1,  2,  0,  0,  0,  0,  0) | PIC: F;
                 ***********************                                          ********

       where the *-indicated part is optional and may be left out. However, in case it is
       present, it is parsed for syntax as well. */

    string str, strCoordName;
    char c;
    vector<int64_t> dummy, GLSM;
    vector<int32_t> GLSMranged;

    bool bPIC = false;
    bool bGLSM = false;

    if (!tok.GetNextToken().GetWord(str)) return false;
    if (str.compare("vertex") != 0) return false;

    // Read in the vertex name and try to insert
    if (!tok.GetNextToken().GetWord(strCoordName)) return false;
    pair< map< string, vector<int32_t> >::iterator, bool> newcoord = mCoords.insert(pair<string, vector<int32_t> >(strCoordName, GLSMranged));
    if (newcoord.second == false)
    {
        ERR_OUT("Double vertex name '" << str << "' used.");
        return false;
    }

    if (!tok.GetNextToken().GetSymbol(c)) return false;
    if (c == '=')
    {
        // Read and discard the vertex data
        if (!tok.GetIntegerList(dummy)) return false;

        if (!tok.GetNextToken().GetSymbol(c)) return false;
    }
    while (c == '|')
    {
        if (!tok.GetNextToken().GetWord(str)) return false;
        if (str.compare("PIC") == 0)
        {
            if (bPIC)
            {
                ERR_OUT("Double usage of the identifier PIC.");
                return false;
            }
            bPIC = true;

            // Read and discard Picard generator
            if (!tok.GetNextToken().GetSymbol(c)) return false;
            if (c != ':') return false;

            if (!tok.GetNextToken().GetWord(str)) return false;

            if (!tok.GetNextToken().GetSymbol(c)) return false;
        }
        else if (str.compare("GLSM") == 0)
        {
            if (bGLSM) return false; // Only once!
            bGLSM = true;

            // Read in the GLSM charges
            if (!tok.GetNextToken().GetSymbol(c)) return false;
            if (c != ':') return false;

            GLSM.clear();
            if (!tok.GetIntegerList(GLSM)) return false;
            size_t nCurGLSMch = GLSM.size();
            if (nCurGLSMch < 1)
            {
                ERR_OUT("Invalid number of GLSM charges supplied.");
                return false;
            }

            // Check that we have a equal number of GLSM charges from all vertices
            if (numGLSMch < 1)
                numGLSMch = nCurGLSMch;
            if (numGLSMch != nCurGLSMch)
            {
                ERR_OUT("Unequal number of GLSM charges.");
                return false;
            }

            if (!tok.GetNextToken().GetSymbol(c)) return false;
        }
        else
            return false;
    }

    if (!bGLSM)
    {
        ERR_OUT("No GLSM data was supplied.");
        return false;
    }

    if (c != ';') return false;

    // We have to check if the GLSM charges are in range and convert to 32-bit vars
    GLSMranged.resize(numGLSMch);
    for (size_t i=0; i<numGLSMch; i++)
    {
        if ((GLSM[i] > MAX_GLSM_RANGE) || (GLSM[i] < -MAX_GLSM_RANGE))
        {
            ERR_OUT("Target divisor charge value " << GLSM[i] << " is out of allowed range.");
            return false;
        }
        GLSMranged[i] = (int32_t) GLSM[i];
    }

    // If we arrive here, then the command line was sane and we can copy the retrieved GLSM data
    newcoord.first->second.clear();
    newcoord.first->second.assign(GLSMranged.begin(), GLSMranged.end());

    return true;
}


bool CInputFile::ParseSRidealCommand(CTokenizer &tok)
{
    /* This function parses a 'srideal' command, for which the maximal prototype syntax line is

       srideal [u1*u2,  u3*u4];

       where all the product variables have to be present in the list of coordinates/vertices.
       This more or less suggest using the 'srideal' command AFTER the 'vertex' commands. */

    vector<string> vec;
    string str;
    char c;

    if (!tok.GetNextToken().GetWord(str)) return false;
    if (str.compare("srideal") != 0) return false;

    if (!tok.GetNextToken().GetSymbol(c)) return false;
    if (c != '[') return false;

    // Work through the array of Stanley-Reisner generators
    while (tok.GetNextToken().GetWord(str))
    {
        vec.clear();

        // Read in a single SR generator
        while (tok.GetNextToken().GetSymbol(c))
        {
            // Check if the variable is in the list of coordinates
            if (mCoords.find(str) == mCoords.end())
            {
                ERR_OUT("Coordinate/vertex '" << str << "' not specified.");
                return false;
            }
            else
            {
                if (vec.size() >= 63)
                {
                    // In principle, this cannot really happen, as we are only allowing for 63 coordinates
                    // in total, i.e. someone must try something like 'srideal [u*u*u*u*u*...];'
                    ERR_OUT("Only a product of up to 63 coordinates is supported per generator");
                    return false;
                }
                vec.push_back(str);
            }

            if (c == ',') break;    // Ends the current generator
            if (c == ']') break;    // Ends the list of generators
            if (c != '*') return false;

            // Another coordinate following
            if (!tok.GetNextToken().GetWord(str)) return false;
        }

// This check has been removed as it sometimes seems to be sensible to use single generators as well...
/*      if (vec.size() < 2)
		{
			// No sensible toric variety has Stanley-Reisner generators consisting of only a single coordinate
			ERR_OUT("Each Stanley-Reisner ideal generator must be the product of at least 2 coordinates.");
			return false;
		}*/

        vSRideal.push_back(vec);

        if (c == ']') break;
        if (c != ',') return false;
    }

    if (!tok.GetNextToken().GetSymbol(c)) return false;
    if (c != ';') return false;

    return true;
}


bool CInputFile::ParseAmbientCohomCommand(CTokenizer &tok)
{
    /* This function parses a 'ambientcohom' cmd, for which the maximal prototype syntax line is

       ambientcohom O(3,0,2,-3,-4,-5,0,2);

       where the number of supplied GLSM charges must comply with the number of GLSM charges
       specified in the vertex commands. However, in principle this command may come first, which
       then defines the number of GLSM charges for the vertices. */

    string str;
    vector<int64_t> GLSM;
    char c;

    if (!tok.GetNextToken().GetWord(str)) return false;
    if (str.compare("ambientcohom") != 0) return false;

    if (!tok.GetNextToken().GetWord(str)) return false;
    if (str.compare("O") != 0) return false;

    // Read in the GLSM charges
    if (!tok.GetIntegerList(GLSM)) return false;
    size_t nCurGLSMch = GLSM.size();
    if (nCurGLSMch < 1)
    {
        ERR_OUT("Invalid number of GLSM charges supplied.");
        return false;
    }

    // Check if we have the correct number of GLSM charges
    if (numGLSMch < 1)
        numGLSMch = nCurGLSMch;
    if (numGLSMch != nCurGLSMch)
    {
        ERR_OUT("Unequal number of GLSM charges.");
        return false;
    }

    if (!tok.GetNextToken().GetSymbol(c)) return false;
    if (c != ';') return false;

    // Now we have to check the 64-bit GLSM data for range and resort them into
    // a vector of 32-bit variables;
    vector<int32_t> GLSMranged;
    GLSMranged.resize(numGLSMch);
    for (size_t i=0; i<numGLSMch; i++)
    {
        if ((GLSM[i] > MAX_GLSM_RANGE) || (GLSM[i] < -MAX_GLSM_RANGE))
        {
            ERR_OUT("Target divisor charge value " << GLSM[i] << " is out of allowed range.");
            return false;
        }
        GLSMranged[i] = (int32_t) GLSM[i];
    }

    // If we arrive here, everything is sane and we can add the data
    vAmbientCohoms.push_back(GLSMranged);

    return true;
}


bool CInputFile::ParseMonomialFileCommand(CTokenizer &tok)
{
    /* This function parses a 'monomialfile' cmd, for which the two prototype syntax lines are

       monomialfile "filename";
       monomialfile off;

       which allows to change the intermediate monomial filename or turn the monomial file usage
       off alltogether. */

    string str;
    char c;

    if (!tok.GetNextToken().GetWord(str)) return false;
    if (str.compare("monomialfile") != 0) return false;

    const CToken &curtok = tok.GetNextToken();
    if (curtok.GetString(str))
        CCmdLineArguments::SetMonomialFileName(str);
    else if (curtok.GetWord(str))
    {
        if (str.compare("off") == 0)
            CCmdLineArguments::SetUseMonomFile(false);
        else
        {
            ERR_OUT("Could not recognize option '" << str << "'.");
            return false;
        }
    }
    else
        return false;

    if (!tok.GetNextToken().GetSymbol(c)) return false;
    if (c != ';') return false;

    return true;
}


bool CInputFile::DiscardCommand(CTokenizer &tok)
{
    /* This function parses and ignored any command like input, i.e. it looks for a starting WORD
	   token and ignored everything up to the final symbol:

       anycommand (...ignored...);

       This is used for ignored commands for some of our internally developed programs for ease
	   of usage concerning input files. */

    string str;

    // We expect a word as the first token of any command
    if (!tok.GetNextToken().GetWord(str)) return false;

    // Now search for the command end symbol ';'
    while (!tok.GetCurToken().IsEndToken())
    {
        char c;

        if (tok.GetNextToken().GetSymbol(c))
        {
            if (c == ';')
                return true;
        }
    }

    return false;
}


bool CInputFile::ParseInput(const char *pInputData)
{
    /* This functions takes responsibility of reading in the data, i.e. it starts the tokenizer
       and subsequently calls the relevant command parsing functions defined above. It also
       does some basic input checking in addition to the input checking done in the individual
       command parsing functions. In case of parsing errors it calls the HandleSyntaxError
       function to appropiately show the erroneous position. */

    Clear();

    // Run the input through the tokenizer
    CTokenizer tok;
    if (!tok.TokenizeInputString(pInputData))
        return false;

    string strCommandWord;
    bool bResult = false;
    CToken token = tok.GetCurToken();

    size_t numCmds = 0;
    size_t nVertexCmds = 0, nSRidealCmds = 0, nAmbientCohomCmds = 0, nMonomialFileCmds = 0;

    // Scan through the keywords
    while (token.GetWord(strCommandWord))
    {
        ptrdiff_t iCmdOffset = token.GetInputOffset();
        numCmds++;

        if (strCommandWord.compare("vertex") == 0)
        {
            // 'vertex' command detected
            nVertexCmds++;
            if (nVertexCmds > 63)
            {
                ERR_OUT("Due to internal limitations only up to 63 vertices are supported.");
                bResult = false;
            }
            else
                bResult = ParseVertexCommand(tok);
        }
        else if (strCommandWord.compare("picardgen") == 0)
        {
			bResult = DiscardCommand(tok);
		}
        else if (strCommandWord.compare("srideal") == 0)
        {
            // 'srideal' command detected
            nSRidealCmds++;
            if (nSRidealCmds > 1)
            {
                ERR_OUT("Multiple 'srideal' commands are forbidden.");
                bResult = false;
            }
            else
                bResult = ParseSRidealCommand(tok);
        }
        else if (strCommandWord.compare("ambientcohom") == 0)
        {
            // 'ambientcohom' command detected
            nAmbientCohomCmds++;
            bResult = ParseAmbientCohomCommand(tok);
        }
        else if (strCommandWord.compare("monomialfile") == 0)
        {
            // 'monomialfile' command detected
            nMonomialFileCmds++;
            if (nMonomialFileCmds > 1)
            {
                ERR_OUT("Multiple 'monomialfile' commands are forbidden.");
                bResult = false;
            }
            else
                bResult = ParseMonomialFileCommand(tok);
        }
        else
        {
            ERR_OUT("Unrecognized command '" << strCommandWord << "' in input file");
            bResult = false;
        }

        // Check if there were some problems in the current run
        if (!bResult)
        {
            HandleSyntaxError(numCmds, pInputData + iCmdOffset, tok.GetPrevToken().GetInputOffset() - iCmdOffset);
            return false;
        }

        // The following routine is more of a internal safety check
        token = tok.GetCurToken();
        if (iCmdOffset == token.GetInputOffset())
        {
            // Internal error - we have not advanced one bit  [THIS SHOULD NOT HAPPEN!]
            ERR_OUT("INTERNAL: Not advancing! Terminating at offset " << iCmdOffset << ".");
            return false;
        }
    }

    // A new command always has to start with a WORD, so if we are here we're either at the end of file or wrong input
    if (!token.IsEndToken())
    {
        ERR_OUT("Expecting command.");
        HandleSyntaxError(++numCmds, pInputData + token.GetInputOffset(), 0);
        return false;
    }

    if (nVertexCmds < 1)       { ERR_OUT("No coordinates are specified in the input file."); return false; }
    if (nSRidealCmds < 1)      { ERR_OUT("No Stanley-Reisner ideal generators are specified in the input file."); return false; }
    if (nAmbientCohomCmds < 1) { ERR_OUT("No ambient line bundle charges are specified in the input file."); return false; }

	if (nVertexCmds > CCmdLineArguments::GetMaxVertices())     { ERR_OUT("Maximum allowed number of " << CCmdLineArguments::GetMaxVertices() << " coordinates exceeded."); return false; }
	if (vSRideal.size() > CCmdLineArguments::GetMaxSRgens())   { ERR_OUT("Maximum allowed number of " << CCmdLineArguments::GetMaxSRgens() << " SR generators exceeded."); return false; }
	if (nAmbientCohomCmds > CCmdLineArguments::GetMaxCohoms()) { ERR_OUT("Maximum allowed number of " << CCmdLineArguments::GetMaxCohoms() << " requested cohomologies exceeded."); return false; }

    return true;
}


bool CInputFile::TranslateToInternalData(CInternalData &intdata)
{
    /* This function translates the data obtained from the file to the internal data formal described
       by the CInternalData class. */

    // First clear out the internal data class
    intdata.Clear();

    // Then translate the coordinates
    size_t numCoords = mCoords.size();
    map<string, vector<int32_t> >::const_iterator it = mCoords.begin();
    for (size_t i=0; i<numCoords; i++)
    {
        InternalCoordData icd;

        // Copy the name and GLSM charges
        icd.strName.assign(it->first);
        icd.GLSMch.Clear();
        for (size_t j=0; j<numGLSMch; j++)
            icd.GLSMch.x[j] = it->second[j];

        // Assign the 64-bit variable and compute the complete union
        icd.liVar = 0x1ull << i;
        intdata.liCompleteUnion |= icd.liVar;

        // Update the monomial width for output
        intdata.nMaxMonomWidth += icd.strName.length() + 1;

        intdata.vCoords.push_back(icd);
        it++;
    }
    intdata.nMaxMonomWidth--;

    // Now translate the Stanley-Reisner ideal
    size_t numSRgens = vSRideal.size();
    for (size_t i=0; i<numSRgens; i++)
    {
        uint64_t liCurGen = 0;
        size_t numGens = vSRideal[i].size();
        for (size_t j=0; j<numGens; j++)
        {
            // Now, this is ugly linear scanning, but it won't be a real bottleneck for the programm...
            for (size_t k=0; k<numCoords; k++)
            {
                if (intdata.vCoords[k].strName == vSRideal[i][j])
                {
                    liCurGen |= intdata.vCoords[k].liVar;
                    break;
                }
            }
        }
        intdata.vSRgens.push_back(liCurGen);
    }
    intdata.numGLSMch = numGLSMch;

    // Finally translate the target ambient space bundle GLSM charges
    size_t numAmbCohoms = vAmbientCohoms.size();
    for (size_t i=0; i<numAmbCohoms; i++)
    {
        i32vec64 vec;
        vec.Clear();

        for (size_t j=0; j<numGLSMch; j++)
            vec.x[j] = vAmbientCohoms[i][j];

        intdata.vTargetDivisors.push_back(vec);
    }

    // Make a quick final consistency checks
    if (intdata.GetDimension() < 1)
        return false;

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


CInternalData::CInternalData()
{
    Clear();
}

void CInternalData::Clear()
{
    vCoords.clear();
    vSRgens.clear();
    vTargetDivisors.clear();
    numGLSMch = 0;
    nMaxMonomWidth = 0;
    liCompleteUnion = 0;
}


bool CInternalData::ReadAndParseInputFile(string strFileName, string strAppend)
{
    /* This function simply redirects the Read&Parse task to the CInputFile class while
       providing some elementary output */

	if (CCmdLineArguments::GetVerboseLevel() >= -5)
	{
		if (!strFileName.empty())
		{
			MSG_OUT("Reading in the input file '" << strFileName << "'...");
			MSG_OUT("");
		}
	}

    CInputFile iof;
    if (!iof.ReadAndParseInputFile(strFileName, strAppend, *this))
        return false;

    return true;
}


string CInternalData::Int64ToCoordProduct(uint64_t liProduct, string strSep, string strZeroVal) const
{
    /* This output function takes a 64-bit variable and converts it into human readable form
       using the appropiate coordinate bit masks. */

    string out;

    if (liProduct == 0)
        return strZeroVal;

    // If we have a non-trivial variable, loop throught the coordinates and apply the
    // corresponding bit masks
    size_t numCoords = vCoords.size();
    bool bSep = false;
    for (size_t i=0; i<numCoords; i++)
    {
        if (liProduct & vCoords[i].liVar)
        {
            if (bSep)
                out.append(strSep);
            else
                bSep = true;
            out.append(vCoords[i].strName);
        }
    }

    return out;
}

string CInternalData::Int64ToCoordProductPadded(uint64_t liProduct, string strSep, string strZeroVal) const
{
    /* This output function takes a 64-bit variable and converts it into human readable form
       using the appropiate coordinate bit masks and adds spaces, such that the strings of all
       possible monomials have equal length. */

    string tmp(Int64ToCoordProduct(liProduct, strSep, strZeroVal));
    tmp.resize(nMaxMonomWidth, ' ');
    return tmp;
}

string CInternalData::Int64ToMonomial(uint64_t liProduct) const
{
    /* This output function checks if the output of bits is requested by the command line argument
       and if so adds this output the string obtained from Int64ToCoordProduct */

    if (CCmdLineArguments::GetShowBits())
        return CBits::IntToBinary(liProduct, vCoords.size()) + " = " + Int64ToCoordProduct(liProduct);
    else
        return Int64ToCoordProduct(liProduct);
}

string CInternalData::Int64ToMonomialPadded(uint64_t liProduct) const
{
    /* This output function checks if the output of bits is requested by the command line argument
       and if so adds this output the string obtained from Int64ToCoordProduct. It also applies
       padding to the string to take the variable length into account. */

    if (CCmdLineArguments::GetShowBits())
        return CBits::IntToBinary(liProduct, vCoords.size()) + " = " + Int64ToCoordProductPadded(liProduct);
    else
        return Int64ToCoordProductPadded(liProduct);
}


void CInternalData::GetCanonicalDivisor(i32vec64 &candiv_out) const
{
    /* This functions computes the canonical divisor, whose charges are the negative of the
       sum of all coordinate charges. */

    candiv_out.Clear();
    size_t numCoords = vCoords.size();
    for (size_t i=0; i<numCoords; i++)
    {
        for (size_t j=0; j<numGLSMch; j++)
            candiv_out.x[j] -= vCoords[i].GLSMch.x[j];
    }
}


void CInternalData::PrintInternalData()
{
    /* This output functions prints all the data contained in the in CInternalData class in
       formatted form. */

    char buf[1024];

    // First print a list of all coordinates
    size_t numCoords = vCoords.size();
    MSG_OUT("    The described ambient space is of dimension " << GetDimension() << ".");
    MSG_OUT("    There are " << numCoords << " coordinates, each having " << numGLSMch << " GLSM charges:");
    for (size_t i=0; i<numCoords; i++)
    {
        safe_sprintf(buf, sizeof(buf), "        coord %2d:  %3s  |  ", (int) i+1, vCoords[i].strName.c_str());
        string strLine(buf);

        for (size_t j=0; j<numGLSMch; j++)
        {
            safe_sprintf(buf, sizeof(buf), "%3d ", (int) vCoords[i].GLSMch.x[j]);
            strLine.append(buf);
        }

        if (CCmdLineArguments::GetShowBits())
        {
            strLine.append("    |    ");
            strLine.append(CBits::IntToBinary(vCoords[i].liVar, numCoords));
        }

        MSG_OUT(strLine);
    }

    MSG_OUT("");

    // Next comes a list of all the Stanley-Reisner ideal generators
    size_t numSRgens = vSRgens.size();
    MSG_OUT("    There are " << numSRgens << " generators of the Stanley-Reisner ideal:");
    for (size_t i=0; i<numSRgens; i++)
    {
        safe_sprintf(buf, sizeof(buf), "        SRgen %2d:  %s", (int) i+1, Int64ToCoordProductPadded(vSRgens[i]).c_str());
        string strLine(buf);
        if (CCmdLineArguments::GetShowBits())
        {
            strLine.append("    |    ");
            strLine.append(CBits::IntToBinary(vSRgens[i], numCoords));
        }
        MSG_OUT(strLine);
    }

    MSG_OUT("");

    // And finally we print the list of all requested ambient cohomologies
    size_t numAmbCohoms = vTargetDivisors.size();
    if (numAmbCohoms == 1)
        MSG_OUT("    There is " << numAmbCohoms << " ambient space sheaf cohomology requested:");
    else
        MSG_OUT("    There are " << numAmbCohoms << " ambient space sheaf cohomologies requested:");
    for (size_t i=0; i<numAmbCohoms; i++)
    {
        safe_sprintf(buf, sizeof(buf), "        cohom %2d:  H^i(A; O(", (int) i+1);
        string strLine(buf);

        for (size_t j=0; j<numGLSMch; j++)
        {
            if (j>0)
                strLine.append(",");
            safe_sprintf(buf, sizeof(buf), "%4d", (int) vTargetDivisors[i].x[j]);
            strLine.append(buf);
        }

        MSG_OUT(strLine << " ))");
    }
}

string CInternalData::PrintInternalDataAsMathematicaScriptInput()
{
    /* This function is used for debug output. It sort of translates the input/internal data stored
       in CInternalData to the form required by the legacy Mathematica 7 script. This allows for
       great convenience when checking the results of the C++ implementation to the results of the
       Mathematica script. */

    // First the coordinates
    string out("    GeometryData = {\n      (*Coordinates*){");

    size_t num_coords = vCoords.size();
    for (size_t i=0; i<num_coords; i++)
    {
        if (i>0)
            out.append(",");
        out.append(vCoords[i].strName);
    }

    // Then the Stanley-Reisner ideal generators as a list of lists
    out.append("},\n      (*Stanley-Reisner*){");

    size_t num_srs = vSRgens.size();
    for (size_t i=0; i<num_srs; i++)
    {
        if (i>0)
            out.append(",");
        out.append("{");
        out.append(Int64ToCoordProduct(vSRgens[i], ","));
        out.append("}");
    }

    // And finally the GLSM charges, which represent the projective equivalence relations
    out.append("},\n      (*Equivalence Relations*){");

    char buf[64];
    for (size_t i=0; i<num_coords; i++)
    {
        if (i>0)
            out.append(",");
        out.append("{");
        for (size_t j=0; j<numGLSMch; j++)
        {
            if (j>0)
                out.append(",");
            safe_sprintf(buf, sizeof(buf), "%d", (int) vCoords[i].GLSMch.x[j]);
            out += buf;
        }
        out.append("}");
    }
    out.append("}\n    };");

    return out;
}


////////////////////////////////////////////////////////////////////////////////////////////////////

// In order to save computational time, there is the possibility to save a monomial list with the
// computed secondary/remnant cohomology to a file and read it back in. This functionality is
// encapsulated in the class CMonomialFile.

class CMonomialFile
{
  private:
    static const std::string Header;
    static const char DELIM;
    static const uint64_t Version;
    static const std::string Unique;
    static const std::string Ambiguous;
    static const std::string EndFile;

  public:
    static bool WriteMonomialsFile(const CInternalData &id, std::string strFileName, const CMonomialsList &ml);
    static bool ReadMonomialsFile(const CInternalData &id, std::string strFileName, CMonomialsList &ml);
};

// File format constants
const string   CMonomialFile::Header = "§$ralfsalg/monomial_file/start$§";
const char     CMonomialFile::DELIM = ' ';
const uint64_t CMonomialFile::Version = 10001;              // Increase after format change!
const string   CMonomialFile::Unique = "unique_monoms";
const string   CMonomialFile::Ambiguous = "ambiguous_monoms";
const string   CMonomialFile::EndFile = "§$ralfsalg/monomial_file/end$§";

bool CMonomialFile::WriteMonomialsFile(const CInternalData &id, string strFileName, const CMonomialsList &ml)
{
    /* This private/internal function writes the data from a CMonomialsList class as well as
       some basic information of the current geometry (taken from CInternalData) to a file. */

    // Open the file, if it already exists discard contents
    ofstream ofs(strFileName.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);
    if (!ofs.is_open())
    {
        ERR_OUT("Could not open the monomial file '" << strFileName << "' for writing.");
        return false;
    }

    // Write the header string, file format version and variety dimension
    ofs << Header << DELIM << Version << DELIM;
    ofs << id.GetDimension() << DELIM;

    // Write the number of uniquely contributing monomials followed by a list
    // of the bit-mask, the contributing cohomology group and the secondary/remnant factor
    ofs << Unique << DELIM << ml.unique_monoms.size() << DELIM;
    for (map<uint64_t, UniqueContribData>::const_iterator itu = ml.unique_monoms.begin(); itu != ml.unique_monoms.end(); itu++)
        ofs << itu->first << DELIM << itu->second.Cohom.nGroup << DELIM << itu->second.Cohom.nFactor << DELIM;

    // Write the number of ambiguously contributing monomials followed by a list
    // of the bit-mask, the number of potential contributions, followed by a list
    // of the cohomology group and the secondary/remant factor for each possibility
    ofs << Ambiguous << DELIM << ml.ambiguous_monoms.size() << DELIM;
    for (map<uint64_t, AmbiguousContribData>::const_iterator ita = ml.ambiguous_monoms.begin(); ita != ml.ambiguous_monoms.end(); ita++)
    {
        size_t num_cohoms = ita->second.vCohoms.size();
        ofs << ita->first << DELIM << num_cohoms << DELIM;
        for (size_t i=0; i<num_cohoms; i++)
            ofs << ita->second.vCohoms[i].nGroup << DELIM << ita->second.vCohoms[i].nFactor << DELIM;
    }

    // Write the footer string
    ofs << EndFile;

    ofs.close();

    return true;
}

bool CMonomialFile::ReadMonomialsFile(const CInternalData &id, string strFileName, CMonomialsList &ml)
{
    /* This private/internal function read the data from an intermediate monomial file
       into the CMonomialsList class. Some very basic consistency check is carried out. */

    // Open the file for reading
    ifstream ifs(strFileName.c_str());
    if (!ifs.is_open())
        return false;

    // Read & check header string
    string header;
    ifs >> header;
    if (Header.compare(header) != 0) return false;

    // Read & check file version
    uint64_t version;
    ifs >> version;
    if (version != Version) return false;

    ml.Clear();

    // Read & check variety dimension to internal data
    size_t var_dim;
    ifs >> var_dim;
    if (var_dim != id.GetDimension()) return false;

    // Read & check beginning of uniquely contributing monomials list
    string unique;
    ifs >> unique;
    if (Unique.compare(unique) != 0) return false;

    // Read the entire list of uniquely contributing monomials and store data
    // in the CMonomialsList class
    unsigned int num_uniques;
    ifs >> num_uniques;
    for (unsigned int i=0; i<num_uniques; i++)
    {
        uint64_t monomial;
        UniqueContribData ucd;
        ifs >> monomial >> ucd.Cohom.nGroup >> ucd.Cohom.nFactor;
        pair<map<uint64_t, UniqueContribData>::iterator, bool> ret = ml.unique_monoms.insert(pair<uint64_t, UniqueContribData>(monomial, ucd));
        if (ret.second == false)
        {
            // File already contains this unique monomial - this should never happen
            return false;
        }
    }

    // Read & check the beginning of the ambiguously contributing monomials list
    string ambiguous;
    ifs >> ambiguous;
    if (Ambiguous.compare(ambiguous) != 0) return false;

    // Read the entire list of ambiguously contributing monomials and the subsequent list
    // of all possible contributions - store data into the CMonomialsList class
    unsigned int num_ambiguous;
    ifs >> num_ambiguous;
    for (unsigned int i=0; i<num_ambiguous; i++)
    {
        uint64_t monomial;
        AmbiguousContribData acd;
        unsigned int num_cohoms;

        // Read monomial, the number of all possible contributions and the corresponding list
        ifs >> monomial >> num_cohoms;
        acd.vCohoms.clear();
        for (unsigned int j=0; j<num_cohoms; j++)
        {
            CohomContrib cc;
            ifs >> cc.nGroup >> cc.nFactor;
            acd.vCohoms.push_back(cc);
        }

        pair<map<uint64_t, AmbiguousContribData>::iterator, bool> ret = ml.ambiguous_monoms.insert(pair<uint64_t, AmbiguousContribData>(monomial, acd));
        if (ret.second == false)
        {
            // File already contains this ambiguous monomial - should never happen
            return false;
        }
    }

    ifs.close();

    ml.ClearRationals();

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


void CMonomialsList::Clear()
{
    unique_monoms.clear();
    ambiguous_monoms.clear();
}

void CMonomialsList::ClearRationals()
{
    /* This function clears the number of rational functions for each monomial of both
       the lists of uniquely and ambiguously contributions. This effectively allows to
       reuse the same CMonomialList for different computations within the same geometry. */

    for (map<uint64_t, UniqueContribData>::iterator itu = unique_monoms.begin(); itu != unique_monoms.end(); itu++)
    {
        itu->second.nRationals = 0;
        itu->second.nRationalsDual = 0;
    }

    for (map<uint64_t, AmbiguousContribData>::iterator ita = ambiguous_monoms.begin(); ita != ambiguous_monoms.end(); ita++)
    {
        ita->second.nRationals = 0,
        ita->second.nRationalsDual = 0;
    }
}


bool CMonomialsList::AddUniqueContribution(uint64_t liMonomial, unsigned int nCohomGroup, unsigned int nFactor)
{
    /* Adds an uniquely contributing monomial and the data obtained from the
       secondary/remnant cohomology. */

    // Prepare the new monomial entry
    UniqueContribData tmp;
    tmp.Cohom.nGroup  = nCohomGroup;
    tmp.Cohom.nFactor = nFactor;
    tmp.nRationals = tmp.nRationalsDual = 0;

    // Add & check if the monomial is already present - this should never happen
    pair<map<uint64_t, UniqueContribData>::iterator, bool> ret = unique_monoms.insert(pair<uint64_t, UniqueContribData>(liMonomial, tmp));
    if (ret.second == false)
    {
        ERR_OUT("Supposedly unique monomial " << CBits::IntToBinary(liMonomial) << " is already in the list of uniquely contributing monomials.");
        return false;
    }

    return true;
}

bool CMonomialsList::AddAmbiguousContribution(uint64_t liMonomial, const vector<CohomContrib> &vContributions)
{
    /* Adds an ambiguously contributing monomial and the data obtained from the
       secondary/remnant cohomology, i.e. the list of all potential cohomology groups where
       this monomial might contribute and the corresponding contribution factor. */

    // Prepare the new ambiguous monomial entry
    AmbiguousContribData tmp;
    tmp.vCohoms.assign(vContributions.begin(), vContributions.end());
    tmp.nRationals = tmp.nRationalsDual = 0;

    // Add & check if the monomial is already presen - this should never happen
    pair<map<uint64_t, AmbiguousContribData>::iterator, bool> ret = ambiguous_monoms.insert(pair<uint64_t, AmbiguousContribData>(liMonomial, tmp));
    if (ret.second == false)
    {
        ERR_OUT("Supposedly unique monomial " << CBits::IntToBinary(liMonomial) << " is already in list of ambiguously contributing monomials.");
        return false;
    }

    return true;
}


bool GetIntSequenceStartEnd(const vector<uint64_t> &seq, size_t &min_out, size_t &max_out)
{
    /* This little helper functions takes an arbitrary length vector of integer variables and
       finds the first and last non-zero entry. */

    bool zero_seq = true;
    size_t sequencelen = seq.size();

    // Loop FORWARDS from the beginning to the end of the vector and find first non-zero element
    for (size_t i=0; i<sequencelen; i++)
    {
        if (seq[i] != 0)
        {
            min_out = i;
            zero_seq = false;
            break;
        }
    }

    // If no such element could be found (i.e. a zero sequence), return false
    if (zero_seq)
    {
        min_out = 1;
        max_out = 0;
        return false;
    }

    // Otherwise loop BACKWARDS from the vector end to the first non-zero element and find the
    // last non-zero element
    for (size_t i=sequencelen-1; i>=min_out; i--)
    {
        if (seq[i] != 0)
        {
            max_out = i;
            break;
        }
    }

    return true;
}

typedef struct
{
    uint64_t liMonomial;
    const UniqueContribData *pUCD;
} OutputUniqueSort;

void CMonomialsList::PrintMonomialList(const CInternalData &id, bool bPrintFactors, bool bPrintRationals, bool bShortList) const
{
    /* This output functions prints an easily readable list of all the data contained in
       the monomial list. The output of the secondary/remnant cohomology factors and the
       output of the potentially already computed number of rational functions can be
       turned on and of. Note that output of the rational functions automaticall turns
       on the output of the secondary/remnant cohomology factors. */

    // We need small helper structure in order to sort the uniquely contributing monomials
    // by their respective cohomology group

    if (bPrintRationals)
        bPrintFactors = true;

    // Store the complete union
    uint64_t complete_union = id.GetCompleteUnion();

    // Presort the unique monomials list
    size_t var_dim = id.GetDimension();
    vector< vector<OutputUniqueSort> > unique_sorted;
    unique_sorted.resize(var_dim+1);

    map<uint64_t, UniqueContribData>::const_iterator itu = unique_monoms.begin();
    while (itu != unique_monoms.end())
    {
        OutputUniqueSort tmp;
        tmp.liMonomial = itu->first;
        tmp.pUCD = &itu->second;
        unique_sorted[itu->second.Cohom.nGroup].push_back(tmp);
        itu++;
    }

    // Print out the unique monomials
    if (unique_monoms.size() < 1)
        MSG_OUT("    There are no unique contribution monomials to the " << var_dim << "-dimensional variety.");
    else
    {
        MSG_OUT("    The " << unique_monoms.size() << " unique contribution monomials to the " << var_dim << "-dimensional variety are:");
        for (unsigned int dim = 0; dim <= var_dim; dim++)
        {
            char buf[256];
            size_t curmonoms = unique_sorted[dim].size();
            MSG_OUT("        H^" << dim << "(X) has " << curmonoms << " contributing unique monomials:");
            for (size_t i=0; i<curmonoms; i++)
            {
                buf[0] = 0;
                uint32_t rat = unique_sorted[dim][i].pUCD->nRationals;
                uint32_t drat = unique_sorted[dim][i].pUCD->nRationalsDual;

				// If a short list is requested, only print the monomials with non-zero contribution
				if (bShortList)
				{
					if ((rat == 0) && (drat == 0))
						continue;
				}

				// Add the secondary/remnant cohomology factors if requested
                if (bPrintFactors)
                {
                    uint32_t factor = unique_sorted[dim][i].pUCD->Cohom.nFactor;
                    // Add the number of rational functions if requested
                    if (bPrintRationals)
                        safe_sprintf(buf, sizeof(buf), " factor %3d * %-3d rationals = contribution %3d  |  factor %3d %-3d dual rationals = contribution %d", (int) factor, (int) rat, (int) (factor * rat), (int) factor, (int) drat, (int) (factor * drat));
                    else
                    {
                        if (factor != 1)
                            safe_sprintf(buf, sizeof(buf), "   factor %3d", (int) factor);
                    }
                }
                MSG_OUT("            " << id.Int64ToMonomialPadded(unique_sorted[dim][i].liMonomial).c_str() << buf);
            }
        }
    }

    // Print out the ambiguous monomials
    if (ambiguous_monoms.size() < 1)
        MSG_OUT("    There are no ambiguous contribution monomials to the variety.");
    else
    {
        MSG_OUT("    There are also " << ambiguous_monoms.size() << " ambiguous monomials:");
        for (map<uint64_t, AmbiguousContribData>::const_iterator ita = ambiguous_monoms.begin(); ita != ambiguous_monoms.end(); ita++)
        {
            char buf[256];
            size_t len = 0;

            size_t numcohoms = ita->second.vCohoms.size();
            for (size_t i=0; i<numcohoms; i++)
            {
                string strtmp;

                if (len < 1)
                {
                    safe_sprintf(buf, sizeof(buf), "            %s ", id.Int64ToMonomialPadded(ita->first).c_str());
                    strtmp = buf;
                    len = strtmp.length();
                }
                else
                {
                    strtmp.assign(len, ' ');
                }

                buf[0] = 0;
                // Add the secondary/remnant cohomology factors if requested
                if (bPrintFactors)
                {
                    uint32_t group = ita->second.vCohoms[i].nGroup;
                    uint32_t factor = ita->second.vCohoms[i].nFactor;
                    uint32_t rat = ita->second.nRationals;
                    uint32_t drat = ita->second.nRationalsDual;
                    // Add the number of rational functions if requested
                    if (bPrintRationals)
                        safe_sprintf(buf, sizeof(buf), "factor %3d * %-3d rationals = contribution %3d  |  %-3d dual rationals = contribution %d  -->  h^%-2d", (int) factor, (int) rat, (int) (factor * rat), (int) drat, (int) (factor * drat), group);
                    else
                    {
                        if (factor != 1)
                            safe_sprintf(buf, sizeof(buf), "factor %3d  -->  h^%-2d", (int) factor, (int) group);
                        else
                            safe_sprintf(buf, sizeof(buf), "            -->  h^%-2d", (int) group);
                    }
                }
                MSG_OUT(strtmp << buf << "   (complement: " << id.Int64ToMonomial(~ita->first & complete_union) << ")");
            }
        }
    }
}


bool CMonomialsList::ReadMonomialsFile(const CInternalData &id, std::string strFileName)
{
    /* Wrapper function to access the reading function of the CMonomialFile class. */
    return CMonomialFile::ReadMonomialsFile(id, strFileName, *this);
}

bool CMonomialsList::WriteMonomialsFile(const CInternalData &id, std::string strFileName) const
{
    /* Wrapper function to access the writing function of the CMonomialFile class. */
    return CMonomialFile::WriteMonomialsFile(id, strFileName, *this);
}