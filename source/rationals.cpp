////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  rationals.cpp                                                                                 //
//  =============                                                                                 //
//                                                                                                //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 08.04.2010  File created as rationals.cpp                                             //
//                      Handles the counting of the rational functionsmology which ultimately     //
//                      yields a list of the relevant denominator monomials and their respective  //
//                      multiplicities.                                                           //
//        - 23.02.2010  Added multi-threading support for the computation of secondary sequences. //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "rationals.h"
#include "main.h"
#include "platform.h"


extern "C" {
#include "polylib/polylib64.h"
}


using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////


string CCohomology::GetCohomologyGroupString(const vector<int32_t> &BundleGLSMcharges)
{
    /* This static output function produces the string "H^i(A; O(2,3,2,4))" with the 
       line bundle's divisor charges. */

    char buf[32];
    string strLine = "H^i(A; O(";
    for (size_t j=0; j<BundleGLSMcharges.size(); j++)
    {
        if (j > 0)
            safe_sprintf(buf, sizeof(buf), ",%4d", (int) BundleGLSMcharges[j]);
        else
            safe_sprintf(buf, sizeof(buf), "%4d", (int) BundleGLSMcharges[0]);
        strLine += buf;
    }
    strLine += " ))";
    return strLine;
}

string CCohomology::GetGroupDimensionsString(const std::vector<uint32_t> &CohomologyDimensions)
{
    /* This static output function produces the string "( 3, 2, 4, 5, 2)" of the cohomology 
       group dimensions. */

    char buf[32];
    string strLine = "(";
    for (size_t i=0; i<CohomologyDimensions.size(); i++)
    {
        if (i > 0)
            safe_sprintf(buf, sizeof(buf), ",%4d", (int) CohomologyDimensions[i]);
        else
            safe_sprintf(buf, sizeof(buf), "%4d", (int) CohomologyDimensions[i]);
        strLine += buf;
    }
    strLine += " )";
    return strLine;
}

string CCohomology::GetGroupDimensionsStringNoPadding(const std::vector<uint32_t> &CohomologyDimensions)
{
    /* This static output function produces the string "3,2,4,5,2" of the cohomology 
       group dimensions without any padding as a comma-seperated list. */

    char buf[32];
    string strLine = "";
    for (size_t i=0; i<CohomologyDimensions.size(); i++)
    {
        if (i > 0)
            safe_sprintf(buf, sizeof(buf), ",%d", (int) CohomologyDimensions[i]);
        else
            safe_sprintf(buf, sizeof(buf), "%d", (int) CohomologyDimensions[i]);
        strLine += buf;
    }
    return strLine;
}

string CCohomology::GetCohomologyString(const vector<int32_t> &BundleGLSMcharges, const vector<uint32_t> &CohomologyDimensions)
{
    /* This static output function naively connects the cohomology group string the 
       cohomology group dimension string, yielding "H^i(A; O(2,3,2,4)) = (3,2,4,5,2)" */

    return GetCohomologyGroupString(BundleGLSMcharges) + " = " + GetGroupDimensionsString(CohomologyDimensions);
}

string CCohomology::GetCohomologyString() const
{
    /* This member function checks if the number of possible cohomology group dimensions is
       ambiguous or could not be determined at all and yields the appropiate output. */

    if (IsAmbiguous())
    {
        if (IsNotDetermined())
            return CCohomology::GetCohomologyGroupString(BundleGLSMch) + " could not be determined";
        else
        {
            char buf[64];
            safe_sprintf(buf, sizeof(buf), " has %d ambiguous solutions", (int) CohomologyDims.size());
            return CCohomology::GetCohomologyGroupString(BundleGLSMch) + buf;
        }
    }
    else
        return CCohomology::GetCohomologyString(BundleGLSMch, CohomologyDims[0]);
}


void CCohomology::PrintFullCohomologyMonomialMap(const CInternalData &id, bool ShortList) const
{
    /* This wrapper function outputs the full list of contributing monomials, including the
       secondary/remnant cohomology factors and computed numbers rational functions. Basically,
       it completely shows the individual contributions to the cohomology group dimensions. */

    MSG_OUT("Monomials, cohomology factors and rational functions of " << GetCohomologyString());
    monoms.PrintMonomialList(id, true, true, ShortList);
}


void CCohomology::PrintCohomologies(const vector<CCohomology> &cohomologies)
{
    /* This function outputs a list of all cohomology group dimensions supplied by the
       argument vector. In case of ambiguous results, each possibility is listed. */

    MSG_OUT("Cohomology dimensions:");
    MSG_OUT("======================");

    // Run through the entire argument vector
    size_t nCohoms = cohomologies.size();
    for (size_t k=0; k<nCohoms; k++)
    {
        MSG_OUT("    dim " << cohomologies[k].GetCohomologyString());
        
        // In case of ambiguous results, show all the possible candidate cohomology group dimensions
        if (cohomologies[k].IsAmbiguous())
        {
            const string offset = "        candidate results are:   ";
            string line;
            for (size_t i=0; i<cohomologies[k].CohomologyDims.size(); i++)
            {
                if (i > 0)
                    line = string(offset.length(), ' ');
                else
                    line.assign(offset);
                MSG_OUT(line << " = " << GetGroupDimensionsString(cohomologies[k].CohomologyDims[i]));
            }
        }
        MSG_OUT("");
    }
}


void CCohomology::SummarizeCohomologies(const vector<CCohomology> &cohomologies)
{
	/* This function outputs a list of lists of all cohomology group dimensions supplied by the
	   argument vector. In case of ambiguous results an error message is places, i.e. it only treats
	   non-ambiguous cohomology groups as "valid" output. */
    
	// Run through the entire argument vector
    size_t nCohoms = cohomologies.size();
	bool bAllRight = true;
	string out = "";
    for (size_t k=0; k<nCohoms; k++)
	{
		// Check if this cohomology is ambiguous
        if (cohomologies[k].IsAmbiguous())
        {
			bAllRight = false;
			break;
		}

		// Write a comma-seperated list of the dimensions
		if (k>0)
			out.append(",");
		out.append("{{");
		out.append(GetGroupDimensionsStringNoPadding(cohomologies[k].CohomologyDims[0]));
		out.append("},{");
		size_t nContribDemons = cohomologies[k].ContributingDenoms.size();
		for (size_t r=0; r<nContribDemons; r++)
		{
			if (r > 0)
				out.append(",");
			out.append(cohomologies[k].ContributingDenoms[r]);
		}
		out.append("}}");
	}

	// Produce output
	if (bAllRight)
		MSG_OUT_NOENDL("{True," << out << "}");
	else
		MSG_OUT_NOENDL("{False,\"Ambiguous cohomologies\"");
}


void CCohomology::GetMathematicaCohomologiesList(const vector<CCohomology> &cohomologies, string &out)
{
    /* This function translates the results of the computations, i.e. the dimensions of the
       line bundle sheaf cohomology groups into a form easily comparable to the condensed output
       of the legacy Mathematica 7 script. */

    // First put out the 'LinebundleCohomologyOf' commands, which basically correspond to our
    // 'ambientcohom' commands and specify the line bundle divisor's charges.
    out = "    (*Requested cohomology commands:*)\n";
    char buf[64];
    size_t nCohoms = cohomologies.size();
    for (size_t k=0; k<nCohoms; k++)
    {
        out += "    LinebundleCohomologyOf[{";
        size_t num_glsm = cohomologies[k].BundleGLSMch.size();
        for (size_t i=0; i<num_glsm; i++)
        {
            if (i>0)
                out += ",";
            safe_sprintf(buf, sizeof(buf), "%d", (int) cohomologies[k].BundleGLSMch[i]);
            out += buf;
        }
        out += "}];\n";
    }

    // Next a list of lists containing the computed cohomology dimensions 
    out += "    ResultingCohomologies = {";

    for (size_t k=0; k<nCohoms; k++)
    {
        if (k>0)
            out += ",";
        out += "{";

        if (cohomologies[k].IsAmbiguous())
        {
            if (cohomologies[k].IsNotDetermined())
                out += "-1 (*not determined*)";
            else
                out += "-1 (*ambiguous*)";
        }
        else
        {
            for (size_t i=0; i<cohomologies[k].CohomologyDims[0].size(); i++)
            {
                if (i>0)
                    out += ",";
                safe_sprintf(buf, sizeof(buf), "%d", (int) cohomologies[k].CohomologyDims[0][i]);
                out += buf;
            }
        }

        out += "}";
    }
    out += "}";
}


////////////////////////////////////////////////////////////////////////////////////////////////////


void StringRemoveSpaces(string &stringIn)
{
	size_t pos = 0;
	bool spacesLeft = true;

	while (spacesLeft)
	{
		pos = stringIn.find(" ");
		if(pos != string::npos)
			stringIn.erase(pos, 1);
		else
			spacesLeft = false;
	}
}

void CRationals::PrintPolyLibConstraintMatrix(const CInternalData &id, void *m)
{
    /* This function prints the equalities and inequalities encoded in a PolyLib matrix structure, 
       which is used in the debug output at the highest verbose level. */

    Matrix *Mat = (Matrix *) m;

    char buf[128];
    string tmp, cureqn, mathematica;
    MSG_OUT("    Condition matrix has " << Mat->NbRows << " rows x " << Mat->NbColumns << " columns");
    if ((Mat->NbRows > 0) && (Mat->NbColumns > 0))
    {
        // Run through all rows of the matrix
        for (unsigned int i=0; i<Mat->NbRows; i++)
        {
            safe_sprintf(buf, sizeof(buf), "        line %2d:   ", (int) i+1);
            tmp = buf;

			// Clear the current equation buffer
			cureqn.clear();
			if (i>0)
				mathematica += ",";

			// Check if we have an equality or inequality
			if (Mat->p[i][0] == 0)
            {
                // equality
                tmp += "eq:    ";
                for (unsigned int j=1; j<Mat->NbColumns-1; j++)
                {
                    if (Mat->p[i][j] == 0)
                        safe_sprintf(buf, sizeof(buf), "          ");
                    else
                    {
                        if (Mat->p[i][j] == 1)
                            safe_sprintf(buf, sizeof(buf), "+ %-7s ", id.GetInternalCoordData()[j-1].strName.c_str());
                        else if (Mat->p[i][j] == -1)
                            safe_sprintf(buf, sizeof(buf), "- %-7s ", id.GetInternalCoordData()[j-1].strName.c_str());
                        else if (Mat->p[i][j] < 0)
                            safe_sprintf(buf, sizeof(buf), "- %3ld*%-3s ", (long int) -Mat->p[i][j], id.GetInternalCoordData()[j-1].strName.c_str());
                        else
                            safe_sprintf(buf, sizeof(buf), "+ %3ld*%-3s ", (long int) Mat->p[i][j], id.GetInternalCoordData()[j-1].strName.c_str());
                    }
                    cureqn += buf;
                }

                safe_sprintf(buf, sizeof(buf), " == %3ld", -(long int) Mat->p[i][Mat->NbColumns-1]);
                cureqn += buf;
            }
            else
            {
                // inquality
                tmp += "ineq:  ";
                for (unsigned int j=1; j<Mat->NbColumns-1; j++)
                {
                    if (Mat->p[i][j] == 0)
                        safe_sprintf(buf, sizeof(buf), "          ");
                    else
                    {
                        if (Mat->p[i][j] == 1)
                            safe_sprintf(buf, sizeof(buf), "+ %-7s ", id.GetInternalCoordData()[j-1].strName.c_str());
                        else if (Mat->p[i][j] == -1)
                            safe_sprintf(buf, sizeof(buf), "- %-7s ", id.GetInternalCoordData()[j-1].strName.c_str());
                        else if (Mat->p[i][j] < 0)
                            safe_sprintf(buf, sizeof(buf), "- %3ld*%-3s ", (long int) -Mat->p[i][j], id.GetInternalCoordData()[j-1].strName.c_str());
                        else
                            safe_sprintf(buf, sizeof(buf), "+ %3ld*%-3s ", (long int) Mat->p[i][j], id.GetInternalCoordData()[j-1].strName.c_str());
                    }
                    cureqn += buf;
                }

                safe_sprintf(buf, sizeof(buf), " >= %3ld", -(long int) Mat->p[i][Mat->NbColumns-1]);
                cureqn += buf;
            }
			
			tmp += cureqn;

			StringRemoveSpaces(cureqn);
			mathematica += cureqn;

            MSG_OUT(tmp);
        }

		tmp = "Reduce[{";
		tmp += mathematica;
		tmp += "},Integers]";
		MSG_OUT("    Mathematica cmd:  " << tmp << "");
    }
}


struct CountRationalsData {
	const CInternalData &id;
	//const vector<InternalCoordData> &coords;
	const i32vec64 &TargetDivisor;
	uint32_t &out;
	uint64_t monomial;
	//size_t numCoords;
	//size_t numGLSMch;
	bool bReturn;
	size_t worker_id;

	CountRationalsData(const CInternalData &init_InternalData,
		uint64_t init_monomial,
		const i32vec64 &init_TargetDivisor,
		uint32_t &init_out,
		size_t init_numCoords,
		size_t init_numGLSMch,
		size_t init_worker_id) :
			id(init_InternalData),
			TargetDivisor(init_TargetDivisor),
			out(init_out)
	{
		monomial = init_monomial;
		//numCoords = init_numCoords;
		//numGLSMch = init_numGLSMch;
		worker_id = init_worker_id;
		bReturn = false; 
	}
};


void CountRationalFunctionsWorker(void *p_dat)
{
    /* This is the main function for the counting of rational functions. If practically serves
       as the wrapper from our data format to the matrix input format required by the PolyLib. */

	CountRationalsData *crd = (CountRationalsData *) p_dat;
	if (!crd) { ERR_OUT("Internet pointer error in monomial counting worker."); return; }

	const vector<InternalCoordData> &coords = crd->id.GetInternalCoordData();
	size_t numCoords = crd->id.GetNumCoordinates();
	size_t numGLSMch = crd->id.GetNumGLSMch();

    // First we need to create a matrix in the PolyLib format
    Matrix *P = Matrix_Alloc((unsigned int) (numCoords + numGLSMch), (unsigned int) (numCoords + 2));
	if (!P) { crd->bReturn = false; ERR_OUT("Internal allocation error in monomial counting worker."); return; }
    const size_t lastcol = numCoords + 2 - 1;

    i32vec64 constraint_vals;
    constraint_vals.Clear();
    for (size_t j=0; j<numGLSMch; j++)
        constraint_vals.x[j] = crd->TargetDivisor.x[j];

    // Now insert the appropiate values and signs into the matrix elements
    for (size_t i=0; i<numCoords; i++)
    {
        if (crd->monomial & coords[i].liVar)
        {
            // Coordinate is in denominator
            for (size_t j=0; j<numGLSMch; j++)
            {
                value_set_si(P->p[j][i+1], - coords[i].GLSMch.x[j]);
                constraint_vals.x[j] += coords[i].GLSMch.x[j];
            }
        }
        else
        {
            // Coordinate is in nominator
            for (size_t j=0; j<numGLSMch; j++)
                value_set_si(P->p[j][i+1], coords[i].GLSMch.x[j]);
        }
    }

    // The first column of the matrix selects the type of the condition
    for (size_t j=0; j<numGLSMch; j++)
    {
        // Coordinates are equalities, so zero the first column
        value_set_si(P->p[j][0], 0);
        // Put the final target charges in last column
        value_set_si(P->p[j][lastcol], -constraint_vals.x[j]);
    }

    // And finally we require the positivity condition on all variables
    for (size_t i=0; i<numCoords; i++)
    {
        value_set_si(P->p[numGLSMch+i][0], 1);
        value_set_si(P->p[numGLSMch+i][i+1], 1);
    }

    // The "universe matrix" is required by PolyLib (no idea what for....)
    Matrix *C = Matrix_Alloc(1, 2);
    value_set_si(C->p[0][0],1); value_set_si(C->p[0][1],1);

	/*  THIS VERBOSE OUTPUT WAS SHIFTED BEHIND THE COMPUTATION
    if (ccmdlinearguments::getverboselevel() >= 5)
    {
        string tmp = "verbose level 6: rationals counting condition matrix for " + id.int64tomonomial(monomial) + ":";
        msg_out(tmp);

        string tmp2(tmp.length(), '-');
        msg_out(tmp2);

        printpolylibconstraintmatrix(id, p);
    }  
	*/

    // Compute a polynomial approximation of the Ehrhart polynomial and evaluate
    Matrix *Validity_Lattice = 0;
    Enumeration *e = Ehrhart_Quick_Apx(P, C, &Validity_Lattice, 0);

    // Now we need to retrieve the actually interesting data, which is at the moment completely
    // idiotically 'solved'. Why is there no function to simply extract the constant part of an e_value?
    // Furtunately, this does not present a bottleneck to the computational speed, but it's ugly like hell...
	long long int val = 0;
	if (e)
	{
		// Brrr.... temporary files for data retrieval... Oh - and we can't use tmpfile() because
		// Microsoft in their infinite wisdom "broke" the command in Windows Vista / 7, as the default
		// destination is only accessible in Administrator-mode...
		char tmpfilename[64];
		safe_sprintf(tmpfilename, sizeof(tmpfilename), "__x_TMP_polylib_file_%dy.tmp", crd->worker_id);
		//const char *tmpfilename = "__x_TMP_polylib_file_y.tmp";
		FILE *points_val = fopen(tmpfilename, "w+");
		if (!points_val)
		{
			ERR_OUT("Could not open temporary file for computation.");
			crd->bReturn = false; return;
		}

		for (Enumeration *en=e; en; en = en->next)
			print_evalue(points_val, &en->EP, NULL);
		rewind(points_val);

		char buf[16];
		int readres = fscanf(points_val, "%c %lld", buf, &val);
		fclose(points_val);
		remove(tmpfilename);
		if (readres < 2)
		{
			ERR_OUT("INTERNAL: Failed to scan the PolyLib result.");
			crd->bReturn = false; return;
		}
	}
	else
	{
		WARN_OUT("Counting of the rationoms errorneous - is your input geometry valid?");

		// During a normal program run, we have the warning output, but in integration mode we treat this as a serious error
		if (CCmdLineArguments::GetVerboseLevel() < -5)
			crd->bReturn = false; return;
	}

    // The variable 'val' now holds our precious number of rational functions
    // A result -1 indicates some serious error in the constraints and therefore in
    // the counting of the rational functions / points in the corresponding polyhedron
    if (val < 0)
    {
        ERR_OUT("INTERNAL: Invalid number of rational functions computed.");
        crd->bReturn = false; return;
    }

    // Debug Output
/*    if (CCmdLineArguments::GetVerboseLevel() >= 5)
    {
		string tmp;
		bool bPrint = true;

		if (val != 0)
			tmp = "verbose level 5: rationals counting condition matrix for " + crd->id.Int64ToMonomial(monomial) + ":";
		else if (CCmdLineArguments::GetVerboseLevel() >= 6)
			tmp = "verbose level 6: rationals counting condition matrix for " + id.Int64ToMonomial(monomial) + ":";
		else
			bPrint = false;

		if (bPrint)
		{
			MSG_OUT(tmp);

			string tmp2(tmp.length(), '-');
			MSG_OUT(tmp2);

			PrintPolyLibConstraintMatrix(id, P);

			MSG_OUT("    Yields value: " << val);
			MSG_OUT("");
		}
    }*/

    // Clean up stuff
    Matrix_Free(C);
    Matrix_Free(P);
    while (e)
    {
        free_evalue_refs(&(e->EP));
        Polyhedron_Free(e->ValidityDomain);
        Enumeration *en = e->next;
        free(e);
        e = en;
    }

    // Store the computed value
	crd->out = (uint32_t) val;

	crd->bReturn = true;
}

std::vector< std::pair<tthread::thread *, void *> > CRationals::WorkersList;

//#define USE_MULTITHREADED_MONOMIAL_COUNTING
// It appears that the PolyLib is NOT thread-safe (causing a weird kind of heap corruption), therefore for the moment we are stuck with a single-threaded evaluation of the rational counting... :-(

#ifdef USE_MULTITHREADED_MONOMIAL_COUNTING

bool CRationals::CountRationalFctns(const CInternalData &id, uint64_t monomial, const i32vec64 &TargetDivisor, uint32_t &out)
{
	CountRationalsData *crd = new CountRationalsData(id, monomial, TargetDivisor, out, id.GetInternalCoordData().size(), id.GetNumGLSMch(), WorkersList.size());
	if (!crd)
		return false; 

	tthread::thread *t = new tthread::thread(CountRationalFunctionsWorker, (void*) crd);
	WorkersList.push_back(make_pair(t, crd));

	return true;
}

bool CRationals::WaitForAllWorkersFinish()
{
	while (WorkersList.size() > 0)
	{
		WorkersList[0].first->join();
		if (!(((CountRationalsData *) WorkersList[0].second)->bReturn))
			return false;
		
		delete WorkersList[0].second;  // free the thread data memory
		delete WorkersList[0].first;   // delete the thread
		WorkersList.erase(WorkersList.begin());
	}
	return true;
}

#else

bool CRationals::CountRationalFctns(const CInternalData &id, uint64_t monomial, const i32vec64 &TargetDivisor, uint32_t &out)
{
	CountRationalsData crd = CountRationalsData(id, monomial, TargetDivisor, out, id.GetInternalCoordData().size(), id.GetNumGLSMch(), 0);

	CountRationalFunctionsWorker(&crd);

	return crd.bReturn;
}

bool CRationals::WaitForAllWorkersFinish()
{
	return true;
}

#endif

void SplitAddToCohomVector(vector<ui32vec64> &inout, const vector<CohomContrib> &split, uint32_t rational, size_t dim)
{
    /* This internal helper functions carries out the branching required for the ultimate attempt to
       resolve the ambiguous contribution. It takes a vector of possible cohomology group dimensions and
       a vector of fixes for a particular contribution. Then all those different candidates are applied to
       all possibilities, which lets the list of all possibilities grow pretty fast. Currently, this is a
       somewhat of a big problem, as it represents the only real "memory hole" in the entire program. */

    size_t inoutlen = inout.size();
    size_t splitlen = split.size();

    ui32vec64 newcohom;
    newcohom.Clear();

    for (size_t i=0; i<inoutlen; i++)
    {
        // Apply the 0th split afterwards to the original elements in inout,
        for (size_t j=1; j<splitlen; j++)
        {
            for (size_t k=0; k<=dim; k++)
                newcohom.x[k] = inout[i].x[k];

            newcohom.x[split[j].nGroup] += split[j].nFactor * rational;

            inout.push_back(newcohom);
        }

        // Finally modify the original in vector
        inout[i].x[split[0].nGroup] += split[0].nFactor * rational;
    }
}

bool CRationals::ComputeCohomology(const CInternalData &id, CMonomialsList &ml, const i32vec64 &TargetDivisor, vector<uint32_t> &out_cohomology, vector<string> &out_contribcohoms)
{
    /* Now, this is the big fat mama of all functions. It takes the monomial list and the line bundle's divisor
       charges as input and computes the corresponding dimensions of the cohomology. Now, we have to keep in mind
       that there still may be ambiguous contributions left over. First we simple count the number of rational
       functions for all those ambiguous cases. If all goes well, none of these troublesome monomials actually
       contributes for real - we can then compute the uniquely determined contributions and are all done. This is
       the nice case. The second attempt is to compute the ambiguous contributions for the Serre-dual divisor and
       hope that all of these are actually zero - then we can simply take the dimensions from the Serre-dual and
       are done. In the worst case, neither the "normal" nor "Serre-dual" ambiguities vanish all at the same time.
       Then we have to take the ugly branch approach, where we compute the set of ALL possible cohomology group
       configurations both for the "normal" and Serre-dual case and check if there are compatible pairs. This may
       ultimately lead to the ambiguous or undeterminant cases. Luckily, this does not happen very often.... */

    size_t dim = id.GetDimension();

    out_cohomology.clear();
    out_cohomology.resize(dim+1);
	out_contribcohoms.clear();

    // Now compute the ambiguous contributions - if any is non-zero we have to take more involved steps
    bool bIsUnique = true;
    for (map<uint64_t, AmbiguousContribData>::iterator it = ml.ambiguous_monoms.begin(); it != ml.ambiguous_monoms.end(); it++)
        if (!CountRationalFctns(id, it->first, TargetDivisor, it->second.nRationals)) return false;
	if (!WaitForAllWorkersFinish())	return false;

    for (map<uint64_t, AmbiguousContribData>::iterator it = ml.ambiguous_monoms.begin(); it != ml.ambiguous_monoms.end(); it++)
    {
        if (it->second.nRationals != 0)
            bIsUnique = false;
    }

    if (bIsUnique)
    {
        // So, there a no contributing non-ambiguous factors - Whohoo!
        // Then compute the unique ones and output
        for (map<uint64_t, UniqueContribData>::iterator it = ml.unique_monoms.begin(); it != ml.unique_monoms.end(); it++)
            if (!CountRationalFctns(id, it->first, TargetDivisor, it->second.nRationals)) return false;
		if (!WaitForAllWorkersFinish()) return false;

        for (map<uint64_t, UniqueContribData>::iterator it = ml.unique_monoms.begin(); it != ml.unique_monoms.end(); it++)
		{
            out_cohomology[it->second.Cohom.nGroup] += it->second.Cohom.nFactor * it->second.nRationals;
			if (it->second.nRationals != 0)
			{
				char buf[32];
				string tmp = "{";
				safe_sprintf(buf, sizeof(buf), "%d,", it->second.Cohom.nGroup);
				tmp += buf;
				safe_sprintf(buf, sizeof(buf), "%d*", it->second.Cohom.nFactor);
				tmp += buf;
				tmp += id.Int64ToCoordProduct(it->first);
				tmp += "}";
				out_contribcohoms.push_back(tmp);
			}
		}

        return true;
    }

    /////////////////////////////////////////////////////////

    // If we have a non-unique situation, compute the rationals for the Serre dual target divisor
    i32vec64 vDualCharges;
    id.GetCanonicalDivisor(vDualCharges);
    size_t numGLSMch = id.GetNumGLSMch();
    for (size_t i=0; i<numGLSMch; i++)
        vDualCharges.x[i] -= TargetDivisor.x[i];

    for (map<uint64_t, UniqueContribData>::iterator it = ml.unique_monoms.begin(); it != ml.unique_monoms.end(); it++)
		if (!CountRationalFctns(id, it->first, vDualCharges, it->second.nRationalsDual)) return false;
	if (!WaitForAllWorkersFinish()) return false;

    bool bDualIsUnique = true;
    for (map<uint64_t, AmbiguousContribData>::iterator it = ml.ambiguous_monoms.begin(); it != ml.ambiguous_monoms.end(); it++)
		if (!CountRationalFctns(id, it->first, vDualCharges, it->second.nRationalsDual)) return false;
	if (!WaitForAllWorkersFinish())	return false;

    for (map<uint64_t, AmbiguousContribData>::iterator it = ml.ambiguous_monoms.begin(); it != ml.ambiguous_monoms.end(); it++)
    {
        if (it->second.nRationalsDual != 0)
            bDualIsUnique = false;
    }

    if (bDualIsUnique)
    {
        // So, there a no contributing non-ambiguous factors in the Serre dual - (little Whohoo!)
        // Compute the unique cases of the Serre-dual, and translate to the "normal" case for output
        for (map<uint64_t, UniqueContribData>::iterator it = ml.unique_monoms.begin(); it != ml.unique_monoms.end(); it++)
		{
            out_cohomology[dim-it->second.Cohom.nGroup] += it->second.Cohom.nFactor * it->second.nRationalsDual;
			if (it->second.nRationalsDual != 0)
			{
				char buf[32];
				safe_sprintf(buf, sizeof(buf), "%d*", it->second.Cohom.nFactor);
				string tmp = buf;
				tmp += id.Int64ToCoordProduct(it->first);
				out_contribcohoms.push_back(tmp);
			}
		}

        return true;
    }

    /////////////////////////////////////////////////////////

    // At this point, it'll get complicated - BOTH the "normal" and "dual" cases are NOT unique...
    // First compute the needed unique normal rationals (those are still left open, when we reach here)
    for (map<uint64_t, UniqueContribData>::iterator it = ml.unique_monoms.begin(); it != ml.unique_monoms.end(); it++)
		if (!CountRationalFctns(id, it->first, TargetDivisor, it->second.nRationals)) return false;
	if (!WaitForAllWorkersFinish()) return false;

    // Now everything is computed (unique/ambiguous contribution for the "normal" and Serre-dual), so we have to
    // figure out how it all fits together.
    // The basic idea is to compute ALL potential output cohomologies, which arise from the ambiguos contributions
    // both on the "normal" and "dual" side and then consider the "intersection" of those cohomologies. Yes, this is
    // very ugly, but at the moment, it seems to be the best we can do...

    // Compute the unique part of all cohomologies, this is basically the starting point for our "ambiguity branching".
    // For speedup, we now use fixed width array of integers... it's somewhat memory wasteful, but considerably faster
    ui32vec64 cohom_dummy, dualcohom_dummy;
    cohom_dummy.Clear();
    dualcohom_dummy.Clear();
    for (map<uint64_t, UniqueContribData>::iterator it = ml.unique_monoms.begin(); it != ml.unique_monoms.end(); it++)
    {
        cohom_dummy.x[it->second.Cohom.nGroup] += it->second.Cohom.nFactor * it->second.nRationals;
        dualcohom_dummy.x[it->second.Cohom.nGroup] += it->second.Cohom.nFactor * it->second.nRationalsDual;
    }
    vector<ui32vec64> cohoms(1, cohom_dummy);
    vector<ui32vec64> dualcohoms(1, dualcohom_dummy);

    // Now loop through the ambiguous part and branch to all possibilities
    // ##############################################
    // ### THIS IS EXTREMELY MEMORY EXPENSIVE!!!! ###  ...and a potential memory hole! IMPROVE!
    // ##############################################
    for (map<uint64_t, AmbiguousContribData>::iterator it = ml.ambiguous_monoms.begin(); it != ml.ambiguous_monoms.end(); it++)
    {
        size_t curCohomsSize = cohoms.size();
        size_t curDualCohomsSize = dualcohoms.size();
        
        // We put the actual branching and check for the looming std::bad_alloc exception, which comes up when we run
        // out of memory. Least we can do in this case is to provide some output of how far we got, but unfortunately,
        // this ends our business...
        try
        {
            if (it->second.nRationals != 0)
			{
                SplitAddToCohomVector(cohoms, it->second.vCohoms, it->second.nRationals, dim);

				// Is this safe? ######################
				//char buf[32];
				//safe_sprintf(buf, "%ud*", it->second.Cohom.nFactor);
				string tmp = "(*Ambiguous*)";
				tmp += id.Int64ToCoordProduct(it->first);
				out_contribcohoms.push_back(tmp);
			}

            curCohomsSize = cohoms.size();

            if (it->second.nRationalsDual != 0)
                SplitAddToCohomVector(dualcohoms, it->second.vCohoms, it->second.nRationalsDual, dim);
        }
        catch(std::bad_alloc &er)
        {
            ERR_OUT("RUNTIME: " << er.what());
            ERR_OUT_PLAIN("Due to the ambiguity branching there are at least " << curCohomsSize << " \"normal\" cohomologies");
            ERR_OUT_PLAIN("and " << curDualCohomsSize << " Serre-dual cohomologies, consuming " << BytesToReadableSize((curCohomsSize + curDualCohomsSize) * sizeof(ui32vec64)) << " memory.");
            ERR_OUT_PLAIN("");
            if (sizeof(void *) < 8)
            {
                ERR_OUT_PLAIN("You could try the same computation again using the 64-bit version of the program.");
                ERR_OUT_PLAIN("");
            }
            return false;
        }
    }

    // OK, the hard part is now done, apparently, we did not run out of memory...
    size_t numcohoms = cohoms.size();
    size_t numdualcohoms = dualcohoms.size();

    if (CCmdLineArguments::GetVerboseLevel() >= 2)
    {
        MSG_OUT("Verbose level 2: Serre dualization information:");
        MSG_OUT("-----------------------------------------------");
        MSG_OUT("    Due to the ambiguous contributions to the requested cohomologies the following braching occured:");
        MSG_OUT("    - \"normal\" configurations:   " << numcohoms << " (" << BytesToReadableSize(numcohoms * sizeof(ui32vec64)) << " memory)");
        MSG_OUT("    - Serre-dual configurations: " << numdualcohoms << " (" << BytesToReadableSize(numdualcohoms * sizeof(ui32vec64)) << " memory)");
        MSG_OUT("");
    }

    // Now find the "Serre-dual intersection" of the computed branching sets
    vector<ui32vec64> stablecohoms;
    for (size_t i=0; i<numcohoms; i++)
    {
        for (size_t j=0; j<numdualcohoms; j++)
        {
            bool bEqual = true;
            for (size_t k=0; k<=dim; k++)
            {
                if (cohoms[i].x[k] != dualcohoms[j].x[dim-k])
                {
                    bEqual = false;
                    break;
                }
            }
            if (bEqual)
                stablecohoms.push_back(cohoms[i]);
        }
    }

    // Clear up the memory of the big branching sets
    cohoms.clear();
    dualcohoms.clear();

    // Now the big decisive moment... could we resolve the ambiguities?
    size_t numstablecohoms = stablecohoms.size();
    if (numstablecohoms != 1)
    {
        // Oh oh.... we were unable to completely resolve the cohomology via Serre
        if (numstablecohoms < 1)
        {
            // No intersection at all? This is certainly bad... and has so far happened not once in testing!
            MSG_OUT("No viable candidate cohomology could be identified via Serre duality.");
            for (size_t i=0; i<=dim; i++)
                out_cohomology[i] = 0;
            // We add a 0 to indicate the number om ambiguous cohomologies
            out_cohomology.push_back((uint32_t) 0);

            return true;
        }
        else
        {
            // Ok, we are not unique, but at least we found some results...
            // Copy the first stable cohomology configuration to the output cohomology vector
            for (size_t i=0; i<=dim; i++)
                out_cohomology[i] = stablecohoms[0].x[i];
            
            // We add a the number of abiguous cohomologies (so n-1 come behind this number)
            out_cohomology.push_back((uint32_t) numstablecohoms);

            // Then we append the other stable cohomologies
            for (size_t k=1; k<numstablecohoms; k++)
            {
                for (size_t i=0; i<=dim; i++)
                    out_cohomology.push_back(stablecohoms[k].x[i]);
            }

            return true;
        }
    }

    // When we reach here, apparently we were able to uniquely resolve the issue
    for (size_t i=0; i<=dim; i++)
        out_cohomology[i] = stablecohoms[0].x[i];

    return true;
}


bool CRationals::ComputeCohomologyAndSerreDual(const CInternalData &id, CMonomialsList &ml, const i32vec64 &TargetDivisor, vector<uint32_t> &out_cohomology, vector<string> &out_contribcohoms)
{
    /* This funcion takes care of the optional Serre-duality check, which is however only carried
       out in case that the results from the "normal" cohomology are unique. Otherwise the Serre-duality
       is already used in the determination of the result cohomologies. */

    if (!ComputeCohomology(id, ml, TargetDivisor, out_cohomology, out_contribcohoms))
        return false;

    // Check if the Serre-duality check is activated...
    if (CCmdLineArguments::GetCheckSerre())
    {
        size_t dim = id.GetDimension();
        if (out_cohomology.size() == dim+1)
        {
            // A unique cohomology was identified, so consider the Serre-dual divisor

            i32vec64 SerreDualDivisor;
            id.GetCanonicalDivisor(SerreDualDivisor);
            size_t numGLSMch = id.GetNumGLSMch();
            for (size_t i=0; i<numGLSMch; i++)
                SerreDualDivisor.x[i] -= TargetDivisor.x[i];

            vector<uint32_t> dual_cohom;
			vector<string> dual_contribcohoms;

            if (!ComputeCohomology(id, ml, SerreDualDivisor, dual_cohom, dual_contribcohoms))
                return false;

            if (dual_cohom.size() == dim+1)
            {
                // If the Serre-dual cohomology could also be determined uniquely, compare
                for (size_t i=0; i<=dim; i++)
                {
                    if (out_cohomology[i] != dual_cohom[dim-i])
                    {
                        // Print a message, if the cohomology does not correspond properly..
                        vector<int32_t> div;
                        div.assign(TargetDivisor.x, TargetDivisor.x + numGLSMch);
                        MSG_OUT("WARNING: The req. cohomology " << CCohomology::GetCohomologyString(div, out_cohomology) << " does not correspond");
                        div.assign(SerreDualDivisor.x, SerreDualDivisor.x + numGLSMch);
                        MSG_OUT("to the Serre dual cohomology " << CCohomology::GetCohomologyString(div, dual_cohom) << ".");
                        return false;
                    }
                }
            }
            else
            {
                vector<int32_t> div;
                div.assign(SerreDualDivisor.x, SerreDualDivisor.x + numGLSMch);
                MSG_OUT("WARNING: The Serre dual cohomology " << CCohomology::GetCohomologyGroupString(div) << " could not");
                MSG_OUT("uniquely be determined. Please check this case manually.");
                return false;
            }
        }
        else
        {
            vector<int32_t> div;
            div.assign(TargetDivisor.x, TargetDivisor.x + id.GetNumGLSMch());
            MSG_OUT("WARNING: The cohomology " << CCohomology::GetCohomologyGroupString(div) << " could not");
            MSG_OUT("uniquely be determined. Please check this case manually.");
            return false;
        }
    }

    return true;
}


bool CRationals::ComputeCohomologies(const CInternalData &id, const CMonomialsList &ml, vector<CCohomology> &out_cohomologies)
{
    /* This function makes a batch run through all the requested line bundle divisors specified
       by the input data. The computed data is then stored in an output vector of CCohomology classes. */

    // Clear the output vector
    out_cohomologies.clear();

    size_t numGLSMch = id.GetNumGLSMch();
    size_t numCohoms = id.GetNumTargetDivisors();

    // Run through all the requested line bundle divisors
    for (size_t i=0; i<numCohoms; i++)
    {
        // Show some progress status output
        char buf[128];
        safe_sprintf(buf, sizeof(buf), "Computing target cohomology %d of %d (%.1f%% done)...", i+1, numCohoms, (double) i*100 / numCohoms);
        CONSOLE_OUT(buf << "       \r");

        // Prepare the input data and compute the cohomology
        const i32vec64 &curtargetdiv = id.GetTargetDivisors()[i];
        
        CCohomology curcohom;
        curcohom.CohomologyDims.clear();
        curcohom.BundleGLSMch.assign(curtargetdiv.x, curtargetdiv.x + numGLSMch);
        curcohom.monoms = ml;
        vector<uint32_t> cohomdims;
		if (!ComputeCohomologyAndSerreDual(id, curcohom.monoms, curtargetdiv, cohomdims, curcohom.ContributingDenoms))
            return false;

        // Sort the data properly into the CCohomology structure
        size_t numentries = cohomdims.size();
        size_t dim = id.GetDimension();
        if (numentries == dim+1)
        {
            // This indicates a non-ambiguous result
            curcohom.CohomologyDims.push_back(cohomdims);
        }
        else if (numentries >= dim+2)
        {
            if (numentries == (cohomdims[dim+1] * (dim+1) + 1))
            {
                // So we have cohomdims[dim+1] results
                vector<uint32_t> tmp(cohomdims.begin(), cohomdims.begin() + dim + 1);
                curcohom.CohomologyDims.push_back(tmp);

                for (size_t j=1; j<cohomdims[dim+1]; j++)
                {
                    tmp.assign(cohomdims.begin() + j*(dim + 1) + 1, cohomdims.begin() + (j+1)*(dim + 1) + 1);
                    curcohom.CohomologyDims.push_back(tmp);
                }
            }
            else if (cohomdims[dim+1] != 0)
            {
                ERR_OUT("Internal error, wrong number of cohomology results.");
                return false;
            }
        }
        else
        {
            ERR_OUT("Internal error, wrong number of cohomology results.");
            return false;
        }

        out_cohomologies.push_back(curcohom);
    }

    CONSOLE_MSG_OUT("Computation of the target cohomology group dimensions complete.");

    return true;
}