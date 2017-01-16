////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  tokenizer.cpp                                                      +----------------------+   //
//  =============                                                      |  generic TOKENIZER   |   //
//                                                                     +----------------------+   //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 04.06.2009  File created as tokenizer.cpp                                             //
//                      Contains a simple generic tokenizer class to read in a string and break   //
//                      it down to its tokens for further parsing.                                //
//        - 17.04.2010  Changed WORD tokens from alpha to alphanums                               //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <stdlib.h>
#include <errno.h>
#include "platform.h"
#include "tokenizer.h"
#include "main.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////


// This macro (yeah, I know... spare me the lecture...) determines which characters we interpret as
// whitespaces, i.e. spaces, tabs, newline and carriage returns
#define IS_WHITESPACE(__x_)  ((__x_ == ' ') || (__x_ == '\t') || (__x_ == '\n') || (__x_ == '\r'))


////////////////////////////////////////////////////////////////////////////////////////////////////


bool operator!=(const CToken &lhs, const CToken &rhs)
{
    /* The operator function handles the comparison of two tokens. We do NOT compare the
       InputOffset, i.e. just the data and type of the tokens has to be equal. */

    if (lhs.tkType != rhs.tkType) return true;

    // Have same type, but do they have same data?
    switch (lhs.tkType)
    {
        case TOKEN_SYMBOL:  return (lhs.symbol != rhs.symbol);
        case TOKEN_INTEGER: return (lhs.integer != rhs.integer);
        case TOKEN_STRING:
        case TOKEN_WORD:    return (lhs.str.compare(rhs.str) != 0);
        case TOKEN_END:
        case TOKEN_ERR:     return false;
    }

    // We should never reach here...
    return false;
}

bool CToken::GetBool(bool &bBool) const
{
    /* As there is no fundamental boolean token, we try to interpret the current token as
       boolean. Therefore, we consider the following values:
        true  -->  WORD:    true (case-inv.)    false  -->  WORD:    false (case-inv.)
                   STRING:  true (case-inv.)           -->  STRING:  false (case-inv.)
                   INTEGER: non-zero value             -->  INTEGER: zero value
       Tokens of any other type cannot be converted to a boolean value and therefore
       return false. */

    switch (tkType)
    {
        case TOKEN_WORD:
        case TOKEN_STRING:
            // We try to recognize either 'true' or 'false' from the string-like tokens
            {
                if (str.length() > 5) return false;
                string strTemp = str;
                transform(strTemp.begin(), strTemp.end(), strTemp.begin(), (int (*) (int)) tolower);
                if ((strTemp.compare("false") == 0) || (strTemp.compare("0") == 0))
                {
                    bBool = false;
                    return true;
                }
                else if ((strTemp.compare("true") == 0) || (strTemp.compare("1") == 0))
                {
                    bBool = true;
                    return true;
                }
                return false;
            }

        case TOKEN_INTEGER:
            // From numerical tokens, we recognize the non-zero values as true
            bBool = (integer != 0);
            return true;

        case TOKEN_SYMBOL:
        case TOKEN_END:
        case TOKEN_ERR:
            // All other tokens cannot be converted to boolean value
            return false;
    }
    return false;
}

string CToken::GetTokenString() const
{
    /* This function converts a token into human readable form, i.e. it prints the token
       type and the token value. */

    char buf[128];
    switch (tkType)
    {
        case TOKEN_SYMBOL:  safe_sprintf(buf, sizeof(buf), "SYMBOL:  %c  (0x%x)", symbol, (int) symbol); break;
        case TOKEN_WORD:    safe_sprintf(buf, sizeof(buf), "WORD:    %s", str.c_str()); break;
        case TOKEN_INTEGER: safe_sprintf(buf, sizeof(buf), "INTEGER: %ld", (long int) integer); break;
        case TOKEN_STRING:  safe_sprintf(buf, sizeof(buf), "STRING:  %s", str.c_str()); break;
        case TOKEN_END:     safe_sprintf(buf, sizeof(buf), "END"); break;
        case TOKEN_ERR:     safe_sprintf(buf, sizeof(buf), "ERR"); break;
        default:            safe_sprintf(buf, sizeof(buf), "--INVALID TOKEN--"); break;
    }
    return buf;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


CTokenizer::CTokenizer()
{
    Clear();
}

void CTokenizer::Clear()
{
    pInputLine = NULL;
    pCurChar = NULL;
    iCurToken = 0;
    vTokens.clear();
}


bool CTokenizer::ReadWord()
{
    /* This internal function tries to read a WORD type token at the current position of the
       input data. A WORD token is specified to be any non-seperated alphanumeric sequence of
       characters not starting with a number. */

    // Determine the length of the alphanumeric sequence
    if (!isalpha(pCurChar[0]))
        return false;
    ptrdiff_t len=1;
    while (isalnum(pCurChar[len]))
        len++;

    // Create a duplicate of the WORD token string
    string strTmp;
    strTmp.assign(pCurChar, len);

    // Create a new token and add to token list
    CToken TmpToken;
    TmpToken.StoreWord(strTmp.c_str(), pCurChar - pInputLine);
    vTokens.push_back(TmpToken);

    pCurChar += len;

    return true;
}

bool CTokenizer::ReadInteger()
{
    /* This internal function tries to read a non-negative INTEGER type token at the current
       position of the input data. Note that the calling function has to take care of a
       potential minus sign in from of this number. If the number does not fit into signed
       64-bit variable, the return value is false. */

    // We determine the length of the numeric sequence
    if (!isdigit(pCurChar[0]))
        return false;
    ptrdiff_t len=1;
    while (isdigit(pCurChar[len]))
        len++;

    // Convert the number
    int64_t iTmp = string_to_int64(pCurChar);
    if (errno == ERANGE)
    {
        // In case the number is larger than 64 bit
        return false;
    }

    // Create a new token and add to the token list
    CToken TmpToken;
    TmpToken.StoreInteger(iTmp, pCurChar - pInputLine);
    vTokens.push_back(TmpToken);

    pCurChar += len;

    return true;
}

bool CTokenizer::ReadString()
{
    /* This internal function tries to read a STRING type token, which is any sequence
       starting and ending with a '"'. Therefore it currently not possible to have a
       '"' character in the string. */

    // Determine the length of the string
    if (pCurChar[0] != '"')
        return false;
    int i=1;
    while ((pCurChar[i] != 0) && (pCurChar[i] != '"'))
        i++;
    if (pCurChar[i] != '"')
        return false;

    // Create a duplicate of the string
    string strTmp;
    strTmp.assign(pCurChar, 1, i-1);

    // Create a new token and add to the token list
    CToken TmpToken;
    TmpToken.StoreString(strTmp.c_str(), pCurChar - pInputLine);
    vTokens.push_back(TmpToken);

    pCurChar += i+1;

    return true;
}

bool CTokenizer::ReadSymbol()
{
    /* This function tries to read in a symbol character. If a '-' character is found it tries
       to read in a subsequent integer and applies the sign. Otherwise the character is simply
       stored as a symbol. Note that this function does NOT check if the character is of
       alphanumeric type, which may allow for a different interpretation. Therefore, this
       function should be the fallback option wenn determining the token type. */

    const char c = pCurChar[0];
    pCurChar++;

    switch (c)
    {
        case '-':
            // If we have a minus sign, look for an integer
            if (ReadInteger())
            {
                // If there is indeed an integer (which is now stored at the last position of
                // the token list) flip the sign
                const size_t numTokens = GetNumberOfTokens();
                int64_t num = 0;
                vTokens[numTokens-1].GetInteger(num);
                vTokens[numTokens-1].StoreInteger(-num, pCurChar - pInputLine - 1);
                break;
            }

        default:
            // If we have no minus sign or cannot find an integer, simply store the character
            // as a SYMBOL token
            CToken TmpToken;
            TmpToken.StoreSymbol(c, pCurChar - pInputLine - 1);
            vTokens.push_back(TmpToken);
    }

    return true;
}

bool CTokenizer::ReadNextToken()
{
    /* This function advances the current character position to the next value which is not
       considered to be a whitespace and then calls the indivual Read*** functions in order
       to properly recognize the token. */

    // First skip all whitespace
    SkipWhitespaces();

    // End of string?
    if (pCurChar[0] == 0)
        return true;

    // Then we read the next token
    if (pCurChar[0] == '"')
        return ReadString();
    if (isdigit(pCurChar[0]))
        return ReadInteger();
    if (isalpha(pCurChar[0]))
        return ReadWord();

    // Note that the ReadSymbol must come last, because everything can be interpreted as an
    // ordinary SYMBOL, which is simply a character.
    return ReadSymbol();
}

void CTokenizer::SkipWhitespaces()
{
    /* Skips whitespace characters (as defined by the macro at the top of the file) until
       either the end of the string of a non-whitespace character is reached. */

    while ((pCurChar[0] != 0) && IS_WHITESPACE(pCurChar[0]))
        pCurChar++;
}

bool CTokenizer::TokenizeInputString(const string &input)
{
    /* This function clears the tokenizer of all prior data and starts the process of
       breaking up the input string into individual tokens. Note that the input string
       is not changed, nor is copy of the original string kept. */

    // Note that the usage of the STL class string ensures that we truly have terminating
    // zero character, i.e. this class is reasonably safe
    pInputLine = input.c_str();
    pCurChar = pInputLine;

    // Read tokens until we hit the end of the string
    while (pCurChar[0] != 0)
    {
        if (!ReadNextToken())
        return false;
    }

    // Add and END token
    CToken EndToken;
    EndToken.SetEndToken(pCurChar - pInputLine);
    vTokens.push_back(EndToken);

    // Set the current token to start
    iCurToken = 0;

    // Clear the "critical" variables for safety
    pInputLine = NULL;
    pCurChar = NULL;

    // OutputTokenList();  // Just for debugging

    return true;
}

void CTokenizer::OutputTokenList() const  // Just for debugging
{
    /* This output function prints a full list of all tokens, which might be useful for
       debugging or writing parsing functions. */

    const size_t numTokens = GetNumberOfTokens();
    MSG_OUT("There are " << numTokens << " tokens in the input string:");
    for (size_t i=0; i<numTokens; i++)
        MSG_OUT("Token " << i << ", Offset " << vTokens[i].GetInputOffset() << ", Type " << vTokens[i].GetTokenString());
}

bool CTokenizer::GetIntegerList(vector<int64_t> &out_list, char cBeginDelim, char cSeperator, char cEndDelim)
{
    /* This function is a semi-parser function, which simplifies the recurring task of reading
       in symbol-seperated integer lists like e.g. comma-seperated bracket-delimited vectors (2,3,1). */

    char c;
    int64_t integer;

    // First we expect a SYMBOL containing the cBeginDelim character
    if (!GetNextToken().GetSymbol(c)) return false;
    if (c != cBeginDelim) return false;

    // The we loop as long as we find and INTEGER followed by the cSeperator character SYMBOL
    while (GetNextToken().GetInteger(integer))
    {
        out_list.push_back(integer);
        if (!GetNextToken().GetSymbol(c)) return false;
        if (c == cEndDelim) break;
        if (c != cSeperator) return false;
    }

    // Finally, there should be a SYMBOL containing the cEndDelim character
    if (c != cEndDelim) return false;

    return true;
}
