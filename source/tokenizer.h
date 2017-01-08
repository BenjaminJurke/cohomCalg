////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  tokenizer.h                                                        +----------------------+   //
//  ===========                                                        |  generic TOKENIZER   |   //
//                                                                     +----------------------+   //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 10.06.2009  File created as tokenizer.h                                               //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef INC_TOKENIZER_H
#define INC_TOKENIZER_H


#include <cstddef>
#include <stdint.h>
#include <vector>
#include <string>


////////////////////////////////////////////////////////////////////////////////////////////////////

// The basic goal of a tokenizer is to break down some input string (i.e. a simple series of characters)
// down into individual tokens (like strings, numbers, special symbols), which can then be easily used 
// for further processing.
//
// The CToken class corresponds to the individual "tokens" to which the tokenizer breaks down the input
// string. 


enum MyTokenType
{
    // Data content
    TOKEN_SYMBOL  = 1,
    TOKEN_WORD    = 2,
    TOKEN_INTEGER = 3,
    TOKEN_STRING  = 4,

    // Control tokens
    TOKEN_END     = 0,
    TOKEN_ERR     = -1
};

class CToken
{
  friend class CTokenizer;

  private:
    MyTokenType tkType;
    ptrdiff_t iInputOffset;

    // Data variables
    char symbol;
    int64_t integer;
    std::string str;

  private:
    inline CToken()
        { tkType = TOKEN_ERR; iInputOffset = 0; symbol = 0; integer = 0; };

    // Data storage/insertion
    inline void StoreSymbol(char cSymbol, ptrdiff_t offset)
         { symbol = cSymbol;   tkType = TOKEN_SYMBOL;  iInputOffset = offset; };
    inline void StoreWord(const char *strWord, ptrdiff_t offset)
         { str = strWord;      tkType = TOKEN_WORD;    iInputOffset = offset; };
    inline void StoreInteger(int64_t iInteger, ptrdiff_t offset)
         { integer = iInteger; tkType = TOKEN_INTEGER; iInputOffset = offset; };
    inline void StoreString(const char *strString, ptrdiff_t offset)
         { str = strString;    tkType = TOKEN_STRING;  iInputOffset = offset; };
    inline void SetEndToken(ptrdiff_t offset)
        { tkType = TOKEN_END; iInputOffset = offset; };

  public:
    // Comparision operator
    friend        bool operator!=(const CToken &lhs, const CToken &rhs);
    friend inline bool operator==(const CToken &lhs, const CToken &rhs) { return (!(lhs != rhs)); };

    // Type and control structures retrival
    inline MyTokenType WhatType() const  
        { return tkType; };
    inline bool      IsEndToken() const     
        { return (WhatType() == TOKEN_END); }
    inline ptrdiff_t GetInputOffset() const 
        { return iInputOffset; }

    // Data retrival
    inline bool GetSymbol(char &cSymbol) const
         { if (tkType == TOKEN_SYMBOL)  { cSymbol = symbol;   return true; } return false; };
    inline bool GetWord(std::string &strWord) const
         { if (tkType == TOKEN_WORD)    { strWord = str;      return true; } return false; };
    inline bool GetInteger(int64_t &iInteger) const
         { if (tkType == TOKEN_INTEGER) { iInteger = integer; return true; } return false; };
    inline bool GetString(std::string &strString) const
         { if (tkType == TOKEN_STRING)  { strString = str;    return true; } return false; };
    bool        GetBool(bool &bBool) const;

    // For debugging purposes
    std::string GetTokenString() const;
};


////////////////////////////////////////////////////////////////////////////////////////////////////


class CTokenizer
{
  private:
    // Data variables
    const char *pInputLine;
    const char *pCurChar;
    std::vector<CToken> vTokens;
    size_t iCurToken;

  private:
    // Internal functions for reading tokens
    bool ReadWord();
    bool ReadInteger();
    bool ReadString();
    bool ReadSymbol();
    bool ReadNextToken();
    void SkipWhitespaces();

  public:
    CTokenizer();
    void Clear();

    // Data retrieval
    inline size_t  GetNumberOfTokens() const          { return vTokens.size(); };
    inline size_t  GetCurTokenIndex() const           { return iCurToken; };
    inline const CToken &GetToken(size_t index) const { return vTokens.at(index); };
    inline const CToken &GetCurToken() const          { return GetToken(iCurToken); };
    inline const CToken &GetNextToken()               { return GetToken(iCurToken++); };      // Note that GetNextToken POST-increments,
    inline const CToken &GetPrevToken()               { return GetToken(--iCurToken); };        // whereas GetPrevToken PRE-decrements!

    bool GetIntegerList(std::vector<int64_t> &out_list, char cBeginDelim = '(', char cSeperator = ',', char cEndDelim = ')');

    // Output functions
    void OutputTokenList() const;

    // Main function to initialize the tokenizer
    bool TokenizeInputString(const std::string &input);
};


#endif
