////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  secondarycohom.cpp                                                                            //
//  ==================                                                                            //
//                                                                                                //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 04.04.2010  File created as secondarycohom.cpp                                        //
//                      Handles the computations of the secondary cohomology which ultimately     //
//                      yields a list of the relevant denominator monomials and their respective  //
//                      multiplicities.                                                           //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iterator>
#include <fstream>
#include <ctime>
#include <cstring>

#include "secondarycohom.h"
#include "main.h"
#include "platform.h"
#include "tinythread.h"


using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////


#define TRAVERSE_SMP_THRESHOLD   8192
#define INSERTION_BUFFERSIZE     (1024*1024)      // Size of the insertion buffer window length
#define MONOMIAL_MAP_MAX_FULL_OUTPUT 250      // Max. number of monomials printed with shorten_output


////////////////////////////////////////////////////////////////////////////////////////////////////


void CSecondaryCohomology::Clear()
{
    insertion_prototype.Clear(0);
    unique_monoms.clear();
    min_c_deg = 0;
    max_c_deg = 0;
    len_c_deg = 0;
    monoms.Clear();
}

CSecondaryCohomology::CSecondaryCohomology()
{
    Clear();
}

string CSecondaryCohomology::Print2ndSeq(const MonomData &mondat, int32_t min_out, int32_t max_out)
{
    /* This output function prints the C-degrees of a secondary/remnant sequence in the
       range C^min_out to C^max_out. The output range values are checked for validity.
       Zero C-degrees in this list are replaced by spaces to make tabular output possible. */

    string strTmp;
    char buf[64];

    // Check the C-degree range
    if (min_out < min_c_deg)
        min_out = min_c_deg;
    if (max_out > max_c_deg)
        max_out = max_c_deg;

    // Loop through the range and print any non-zero C-degree, otherwise spaces
    for (int32_t j=min_out; j<=max_out; j++)
    {
        safe_sprintf(buf, sizeof(buf), "C^%d=%lld  ", (int) j, (long long int) mondat.c[j-min_c_deg]);
        if (mondat.c[j-min_c_deg] != 0)
            strTmp += buf;
        else
            strTmp += string(strlen(buf), ' ');
    }
    return strTmp;
}

string CSecondaryCohomology::Print2ndSeq(const MonomData &mondat)
{
    /* This output functions prints the C-degrees of the secondary/remnant sequence starting
       with the first and ending with the last non-zero C-degree. Note that this is
       NOT suitable for tabular output. */

    size_t seq_start, seq_end;
    if (GetIntSequenceStartEnd(mondat.c, seq_start, seq_end))
        return Print2ndSeq(mondat, (int32_t) seq_start+min_c_deg, (int32_t) seq_end+min_c_deg);
    else
        return "zero sequence";
}

string CSecondaryCohomology::PrintMonomialWith2ndSeq(const CInternalData &id, MonomMap::const_iterator monom, int32_t min_out, int32_t max_out)
{
    /* This output function prints the monomial (and potentially the internal bit-mask), determines
       the length of the non-zero part of the C-degrees and prints the range of C-degrees from
       min_out to max_out. */

    char buf[64];
    size_t seq_start = 0, seq_end = 0;
    safe_sprintf(buf, sizeof(buf), " with %4ld total (len %2d):   ", (long int) monom->second.num, (int) (seq_end - seq_start + 1));

    // Check if we have a zero sequence
    if (GetIntSequenceStartEnd(monom->second.c, seq_start, seq_end))
        return id.Int64ToMonomialPadded(monom->first) + buf + Print2ndSeq(monom->second, min_out, max_out);
    else
        return id.Int64ToMonomialPadded(monom->first) + buf + "zero sequence";
}


void CSecondaryCohomology::PrintMonomMap(const CInternalData &id, bool shorten_output)
{
    /* This output function prints the entire monomial map in tabular form. Only the necessary
       part of the secondary/remnant sequences is printed. Furthermore an automatic shortening
       of excessive large outputs can be activated. */

    // Determine the smalles and highest non-trivial c-degree in the entire list
    int32_t min_deg = len_c_deg;
    int32_t max_deg = 0;
    for (MonomMap::const_iterator set_iter = unique_monoms.begin(); set_iter != unique_monoms.end(); set_iter++)
    {
        int32_t i;

        // Check FORWARDS from the beginning
        for (i = 0; i<min_deg; i++)
        {
            if (set_iter->second.c[i] != 0)
            {
                min_deg = i;
                break;
            }
        }

        // Check BACKWARDS from the end
        for (i = len_c_deg-1; i >= max_deg; i--)
        {
            if (set_iter->second.c[i] != 0)
            {
                max_deg = i+1;
                break;
            }
        }
    }
    min_deg += min_c_deg;
    max_deg += min_c_deg;

    // Now print the monomials
    size_t numelems = unique_monoms.size();
    MonomMap::const_iterator set_iter = unique_monoms.begin();
    MSG_OUT("There are " << numelems << " elements in the monomial map (using " << (len_c_deg * sizeof(int64_t) * numelems + sizeof(int64_t) + sizeof(uint64_t)) << " bytes of memory):");
    for (size_t i=0; i<numelems; i++)
    {
        MSG_OUT("    " << PrintMonomialWith2ndSeq(id, set_iter, min_deg, max_deg));
        set_iter++;

        // Check if the output should be shortened
        if (shorten_output)
        {
            if (numelems > MONOMIAL_MAP_MAX_FULL_OUTPUT)
            {
                if ((i >= 50) && (i < numelems - 50))
                {
                    MSG_OUT("(.....)");
                    i = numelems - 50;
                }
            }
        }
    }
}


bool CSecondaryCohomology::IsSequenceExact(const vector<uint64_t> &seq, size_t &pos_out) const
{
    /* This computational helper function checks a arbitrary length sequence for naive
       exactness. It also computes the naive kernel via the alternating sum method from
       the end and put out a potential repair position in case of non-exactness. */

    size_t seqlen = seq.size();
    if (seqlen > 0)
    {
        int64_t kernel = seq[seqlen-1];
        for (size_t i=0; i<seqlen; i++)
        {
            kernel = seq[seqlen-(i+1)] - kernel;
            if (kernel < 0)
            {
                pos_out = seqlen - (i+1)+1 + min_c_deg;
                return false;
            }
        }
        return true;
    }
    else
    {
        // A zero-length sequence - perhaps throw some runtime error?!
        WARN_OUT("Internal error while checking secondary sequence exactness.");
        return true;
    }
}


/////////////////////////
// The big main functions


bool CSecondaryCohomology::Init2ndSequences(const CInternalData &id)
{
    /* This function initializes the CSecondaryStructure and computes the theoretical
       minimal and maximal C-degree level. It also prepares the corresponding insertion
       prototype for speedup. */

    // Clear the list of unique monomials resulting from unions;
    unique_monoms.clear();

    // Determine the min/max values of the secondary complex
    min_c_deg = 1 - (int32_t) id.GetNumSRgens();
    max_c_deg = (int32_t) id.GetNumCoordinates();
    len_c_deg = max_c_deg - min_c_deg + 1;

    // Prepare the insertion prototype dummy for speedup
    insertion_prototype.Clear(len_c_deg);

    return true;
}


struct FastMonomialData
{
    uint64_t monomial;
    int nka;

  public:
    friend inline bool operator<(const FastMonomialData &lhs, const FastMonomialData &rhs) { return lhs.monomial < rhs.monomial; };
};

void SortAndFlushBuffer(vector<FastMonomialData> &buffer, MonomMap &unique_monoms, size_t max_run, MonomData &insertion_prototype, int min_c_deg)
{
    /* This internal helper function flushes the monomial insertial buffer. As the monomials
       are kept in an associatate list, it helps to presort the momials, such that the number
       of look-ups in the associative index is minimized. */

    if (max_run > 0)
    {
        // Sort the insertion buffer
        sort(buffer.begin(), buffer.begin()+max_run);

        // The first one we have to get by hand
        pair<MonomMap::iterator, bool> ret = unique_monoms.insert(pair<uint64_t, MonomData>(buffer[0].monomial, insertion_prototype));
        ret.first->second.num++;
        ret.first->second.c[buffer[0].nka - min_c_deg]++;

        // Now loop through the buffer
        for (size_t i=1; i<max_run; i++)
        {
            // Check if the monomial is still the same, otherwise get a new one
            if (ret.first->first != buffer[i].monomial)
                ret = unique_monoms.insert(pair<uint64_t, MonomData>(buffer[i].monomial, insertion_prototype));

            // Increase data
            ret.first->second.num++;
            ret.first->second.c[buffer[i].nka - min_c_deg]++;
        }
    }
}


struct TraverseWorkerData
{
	unsigned int worker_id;
	MonomMap unique_monoms;
	vector<FastMonomialData> vBuffer;
	uint64_t window_start, window_end, window_width, buffer_runs;
	MonomData *insertion_prototype;
	int32_t min_c_deg;
	vector<uint64_t> SRunion[4];
	vector<double> *Status;
};

void TraverseWorkerFunction(void *p_dat)
{
	TraverseWorkerData *wd = (TraverseWorkerData *) p_dat;
	FastMonomialData *pbuf;

	// This is the expensive main loop
    for (uint64_t buffer_loop = 0; buffer_loop < wd->buffer_runs; buffer_loop++)
    {
		uint64_t window_start = wd->window_start + buffer_loop * INSERTION_BUFFERSIZE;
        uint64_t window_end = window_start + INSERTION_BUFFERSIZE - 1;
		if (window_end > wd->window_end)
			window_end = wd->window_end;

        pbuf = &wd->vBuffer[0];

        for (uint64_t sr_combi=window_start; sr_combi<=window_end; sr_combi++, pbuf++)
        {
            uint64_t genunion = wd->SRunion[0][sr_combi & 0xffffu] | wd->SRunion[1][(sr_combi >> 16) & 0xffffu] | wd->SRunion[2][(sr_combi >> 32) & 0xffffu] | wd->SRunion[3][(sr_combi >> 48) & 0xffffu];
            pbuf->monomial = genunion;
            pbuf->nka = CBits::CountBits(genunion) - CBits::CountBits(sr_combi);
        }

		SortAndFlushBuffer(wd->vBuffer, wd->unique_monoms, (size_t) (window_end-window_start + 1), *wd->insertion_prototype, wd->min_c_deg);

        // Show some progress status output
		(*wd->Status)[wd->worker_id] = (double) (window_end - wd->window_start) / wd->window_width;
    }
}

struct StatusWorkerData {
	vector<double> *status;
	clock_t start;
	bool Run;
};

void StatusOutputWorker(void *p_dat)
{
	StatusWorkerData *wstat = (StatusWorkerData *) p_dat;
	if (!wstat)
		return;

	while (wstat->Run)
	{
		double dTotal = 0.0;
		for (size_t i=0; i<wstat->status->size(); i++)
			dTotal += (*wstat->status)[i];
		dTotal /= wstat->status->size();
		double dSecElapsed = (double) (clock() - wstat->start) / CLOCKS_PER_SEC;
		double dSecRemaining = (dTotal > 0.0) ? ((1.0 - dTotal) / dTotal) * dSecElapsed : 1.0;
		dTotal *= 100;
		
		char buf[64];
		safe_sprintf(buf, sizeof(buf), "%.2f", dTotal);
		CONSOLE_OUT("  " << buf << "% completed (" << SecondsToTime((clock_t) dSecRemaining) << " remaining)...");

		CONSOLE_OUT("            \r");

		SleepMilliSec(100);
	}

	delete wstat;
}


bool CSecondaryCohomology::TraverseSRpowerset(const CInternalData &id)
{
    /* This is the main powerset traverse function, which loops through the entire range
       of subsets of the Stanley-Reisner generators and computes the corresponding
       monomial unions as well as the secondary/remnant sequences. */

    // As an enourmous speedup, we are using a union buffer analogous to the counting of the bits
    // Since we are dealing with 64 bit variables, we are using 4 bit buffer for 16-bit blocks, such
    // that the buffer size is 4*8*64 KB = 2 MBr. Note that this buffer has to be allocated
    // dynamically, otherwise the MS VS compiler throws a stack overflow in debug mode.


	// Prepare and setup the worker threads
	size_t numWorkers = tthread::thread::hardware_concurrency();
	if (numWorkers == 0)
		numWorkers = 4;

    // Compute the maximal SR bitmask value and the number of buffer runs
    size_t number_of_srgens = id.GetNumSRgens();
    uint64_t sr_powerset_max = (0x1ull << number_of_srgens);
	if (sr_powerset_max < TRAVERSE_SMP_THRESHOLD)
		numWorkers = 1;
	uint64_t sr_powerset_max_per_worker = sr_powerset_max / numWorkers;

    // Allocate the insertion buffer (typically around 12-16 MB)
	TraverseWorkerData *workerdat = new TraverseWorkerData[numWorkers];
	vector<uint64_t> paddedSRvec(id.GetSRgens());
	paddedSRvec.resize(64, 0);
	vector<double> WorkerStatus;
	WorkerStatus.resize(numWorkers, 0.0);
	for (size_t k = 0; k<numWorkers; k++)
	{
		workerdat[k].vBuffer.resize(INSERTION_BUFFERSIZE);
		workerdat[k].window_start = k*sr_powerset_max_per_worker;
		workerdat[k].worker_id = k;
		workerdat[k].window_end = (k+1)*sr_powerset_max_per_worker - 1;
		if (k==numWorkers-1)
			workerdat[k].window_end = sr_powerset_max - 1;
		workerdat[k].window_width = workerdat[k].window_end - workerdat[k].window_start + 1;
		workerdat[k].buffer_runs = workerdat[k].window_width / INSERTION_BUFFERSIZE;
		if (workerdat[k].window_width % INSERTION_BUFFERSIZE != 0)
			workerdat[k].buffer_runs++;
		workerdat[k].insertion_prototype = &insertion_prototype;
		workerdat[k].min_c_deg = min_c_deg;
		workerdat[k].Status = &WorkerStatus;

		// Init the worker buffer
		for (size_t i=0; i<4; i++)
		{
			workerdat[k].SRunion[i].resize(0x1u << 16);
			for (size_t j=0; j<(0x1u << 16); j++)
			{
				uint64_t combi = j << (i*16);
				uint64_t genunion = 0;
				for (size_t k=0; k<16; k++)
				{
					if (combi & (0x1ull << (k + i*16)))
						genunion |= paddedSRvec[k + i*16];
				}
				workerdat[k].SRunion[i][j] = genunion;
			}
		}
		
	}
	paddedSRvec.clear();

	// Start the workers
	vector<tthread::thread *> worker_threads;
	for (size_t k = 0; k<numWorkers; k++)
	{
		tthread::thread *t = new tthread::thread(TraverseWorkerFunction, (void*) &workerdat[k]);
		worker_threads.push_back(t);		
	}
	
	StatusWorkerData *wstat = new StatusWorkerData;
	if (!wstat)
		return false;
	wstat->Run = true;
	wstat->status = &WorkerStatus;
	wstat->start = clock();
	tthread::thread *stat = new tthread::thread(StatusOutputWorker, wstat);

	// Wait for the threads to complete...
	for (size_t k = 0; k<numWorkers; k++)
	{
		worker_threads[k]->join();
		delete worker_threads[k];
	}
	wstat->Run = false;
	

	// Piece together the individually generated data
	unique_monoms.clear();
	for (size_t k = 0; k<numWorkers; k++)
	{
		MonomMap::iterator it = workerdat[k].unique_monoms.begin();

		while (it != workerdat[k].unique_monoms.end())
		{
			pair<MonomMap::iterator, bool> ret = unique_monoms.insert(pair<uint64_t, MonomData>(it->first, it->second));
			if (ret.second==false)
			{
				size_t chain_len = it->second.c.size();
				ret.first->second.num += it->second.num;
				for (size_t r = 0; r < chain_len; r++)
					ret.first->second.c[r] += it->second.c[r];
			}
			it++;
		}
	}

    CONSOLE_OUT(string(75, ' ') << "\r");

    return true;
}


bool CSecondaryCohomology::Compute2ndCohomFromTrivialSequences()
{
    /* After generating the secondary/remnant sequences, the secondary/remnant "cohomology"
       has to be computed from this data. However, lacking a proper description of the mappings
       one has to take a couple of assumptions and employ the Serre-duality to resolve a couple
       of ambiguities. For now, we are taking a two-step scan through the entire range of the
       computed monomials, where in the first one we are simply removing the trivial cases:
       If only a single number is not equal to 0, the sequence is non-exact, but the remnant
       cohomology is obvious. Furthermore, we remove all obviously exact sequence, which is
       quickly checked by computing the alternating sum. */

    MonomMap::iterator cur = unique_monoms.begin();
    while (cur != unique_monoms.end())
    {
        if (cur->second.num == 1)
        {
            // Take care of a sequence containing just a single value 1
            for (size_t i=0; i<len_c_deg; i++)
            {
                if (cur->second.c[i] != 0)
                {
                    if (!monoms.AddUniqueContribution(cur->first, (int32_t) i+min_c_deg, 1))
                        return false;
                    break;
                }
            }
            unique_monoms.erase(cur++);
        }
        else
        {
            // Compute the alternating sum of the secondary/remnant sequence
            int64_t altersum = 0;
            for (size_t i=0; i<len_c_deg; i++)
            {
                if (i % 2)
                    altersum += cur->second.c[i];
                else
                    altersum -= cur->second.c[i];
            }

            if (altersum == 0)
            {
                // Sequence is exact, so discard
                unique_monoms.erase(cur++);
            }
            else
            {
                // Sequence is non-exact & non-trivial, so keep the altersum value AND compute the sequence length
                uint32_t seq_start = len_c_deg;
                for (uint32_t i=0; i<len_c_deg; i++)
                {
                    if (cur->second.c[i] != 0)
                    {
                        seq_start = i;
                        break;
                    }
                }

                uint32_t seq_end = 0;
                for (uint32_t i=len_c_deg-1; i>=0; i--)
                {
                    if (cur->second.c[i] != 0)
                    {
                        seq_end = i;
                        break;
                    }
                }

                uint32_t seq_len = seq_end - seq_start + 1;
                if (seq_len < 1)
                {
                    // This would imply a complete zero sequence!
                    ERR_OUT("Internal error while computing secondary cohomology.");
                    return false;
                }

                // Take care of special length
                switch (seq_len)
                {
                    case 1:
                        // Sequence has a single non-zero value
                        if (!monoms.AddUniqueContribution(cur->first, seq_start+min_c_deg, (uint32_t) cur->second.c[seq_start]))
                            return false;

                        break;

                    default:
                        // Otherwise, store the alternating sum value
                        cur->second.altersum = altersum;
                        cur++;
                }
            }
        }
    }

    return true;
}


bool CSecondaryCohomology::Compute2ndCohomFromSequences(const CInternalData &id)
{
    /* In the second part of the secondary/remnant cohomology computation, we are taking care
       of the non-trivial cases remaining. In principle this procedure is suitable for all
       cases, however, it is somewhat more expensive, therefore the first run to remove the
       "obvious" secondary/remnant sequences. The basic idea is to "repair" the secondary
       sequence in the C-degree range of 0...dim in order to make it exact again. The repair
       value and its possible positions then yield the cohomology of this sequence. On the
       other hand, this is ultimately the origin of the troublesome ambiguities, which can
       only be properly removed by understanding the mappings between the C-degrees. */

    // The temporary "repair" storage
    vector<CohomContrib> fixes;

    // Loop through all remaining monomials in the list
    size_t variety_dim = id.GetDimension();
    MonomMap::iterator cur = unique_monoms.begin();
    while (cur != unique_monoms.end())
    {
        if (cur->second.altersum != 0)
        {
            uint64_t absaltersum = abs(cur->second.altersum);
            fixes.clear();

            // Search for the the C-degrees in range, which is >= the absolute value of the
            // alternating sum of the sequence, which measures the "total failure" exactness
            for (unsigned int j=0; j<=variety_dim; j++)
            {
                if (cur->second.c[j-min_c_deg] >= absaltersum)
                {
                    // Now descrease the value of this position by the absolute value of the
                    // alternating sum and check, if this makes the sequence exact
                    vector<uint64_t> tempseq(cur->second.c);
                    tempseq[j-min_c_deg] -= absaltersum;
                    size_t dummy;
                    if (IsSequenceExact(tempseq, dummy))
                    {
                        // If exactness can be restored this way, the cohomology group at the
                        // respective repair position has the dimension of the absolute value
                        // of the alternating sum - we therefore store this fix as a candidate
                        CohomContrib newfix;
                        newfix.nGroup = j;
                        newfix.nFactor = (uint32_t) absaltersum;
                        fixes.push_back(newfix);
                    }
                }
            }

            // Check if we have obtained at least one "repair" value / secondary cohomology group
            size_t num_fixes = fixes.size();
            if (num_fixes < 1)
            {
                // This should NEVER happen, but is potentially possible, as the
                // repair range is restricted to 0...dim instead of the total sequence
                ERR_OUT("FATAL: Could no properly resolve the exactness of a secondary sequence");
                MSG_OUT(Print2ndSeq(cur->second));
                return false;
            }

            // Now, do we have multi-contributions or not?
            if (num_fixes > 1)
            {
                // Ambiguous contribution we have to take care of lateron
                if (!monoms.AddAmbiguousContribution(cur->first, fixes))
                    return false;
            }
            else
            {
                // Unique contribution... that's fine
                if (!monoms.AddUniqueContribution(cur->first, fixes[0].nGroup, fixes[0].nFactor))
                    return false;
            }
        }
        unique_monoms.erase(cur++);
    }

    // At this point there should be no more monomials left, but check anyway
    if (unique_monoms.size() > 0)
    {
        ERR_OUT("Some monomials could not be properly processed.");
        return false;
    }

    // Clear up allocated buffer space of the monomials (should not be too much...)
    unique_monoms.clear();

    return true;
}


bool CSecondaryCohomology::Perform2ndCohomSerreReduction(const CInternalData &id)
{
    /* Now we have to take care of the ambiguos contributions derive from from restoring the
       exactness of the secondary/remnant sequence. The idea is to check if the complement
       monomial can already be found in the list of unique contributions, because this would
       obviously correspond to a Serre-dual contribution. If so, we can turn the ambiguous
       contribution in an unique one. However, this does not solve every case. In the second
       scan, we check if from the remaining ambiguous contributions an ambiguous complement
       monomial can be found, whose ambiguous contributions might be just right to identify
       a unique one. Unfortunately, this still leaves some cases behind, which we have to deal
       with in the actual computation of the cohomology dimensions. */

    size_t num_ambiguous_before = monoms.ambiguous_monoms.size();
    uint64_t complete_union = id.GetCompleteUnion();
    uint32_t variety_dim = (uint32_t) id.GetDimension();

    // Scan through all ambiguous contributions and check if the complement monomial is unique
    for (map<uint64_t, AmbiguousContribData>::iterator ita = monoms.ambiguous_monoms.begin(); ita != monoms.ambiguous_monoms.end(); )
    {
        uint64_t monomial = ita->first;
        uint64_t complement = (~monomial) & complete_union;

        // First try to find a UNIQUE complement partner
        map<uint64_t, UniqueContribData>::const_iterator uc = monoms.unique_monoms.find(complement);
        if (uc != monoms.unique_monoms.end())
        {
            // Check if this complement partner has a compatible contribution, which would identify
            // a correct contribution of our ambiguous case
            uint32_t group_suggestion = variety_dim - uc->second.Cohom.nGroup;
            size_t numcandidates = ita->second.vCohoms.size();
            bool bFoundCandidate = false;
            size_t candidate = 0;
            for (size_t i=0; i<numcandidates; i++)
            {
                if (ita->second.vCohoms[i].nGroup == group_suggestion)
                {
                    bFoundCandidate = true;
                    candidate = i;
                    break;
                }
            }

            if (bFoundCandidate)
            {
                // So the unique complement monomial exists and identifies the correct contribution
                UniqueContribData newunique;
                newunique.Cohom.nGroup = group_suggestion;
                newunique.Cohom.nFactor = ita->second.vCohoms[candidate].nFactor;
                newunique.nRationals = ita->second.nRationals;
                newunique.nRationalsDual = ita->second.nRationalsDual;

                // Now add the identified unique one
                pair<map<uint64_t, UniqueContribData>::iterator, bool> nu = monoms.unique_monoms.insert(pair<uint64_t, UniqueContribData>(monomial, newunique));
                if (nu.second == false)
                {
                    ERR_OUT("The monomial " << id.Int64ToMonomial(monomial) << " already exists in the unique monomials list!");
                    return false;
                }

                // Remove the ambiguous contribution and start the scan again at the beginning
                monoms.ambiguous_monoms.erase(ita);
                ita = monoms.ambiguous_monoms.begin();
            }
            else
                ita++;
        }
        else
        {
            // If no UNIQUE partner can be found, try to find an ambiguous partner...
            map<uint64_t, AmbiguousContribData>::iterator ita2 = monoms.ambiguous_monoms.find(complement);
            if (ita2 != monoms.ambiguous_monoms.end())
            {
                // OK, we have a partner. Then try to find a "Serre-dual" intersection of the potential contributions
                vector<CohomContrib> intersection;
                size_t cohoms_ita  = ita->second.vCohoms.size();
                size_t cohoms_ita2 = ita2->second.vCohoms.size();
                intersection.clear();
                for (size_t i=0; i<cohoms_ita; i++)
                {
                    for (size_t j=0; j<cohoms_ita2; j++)
                    {
                        if (ita->second.vCohoms[i].nGroup == (variety_dim - ita2->second.vCohoms[j].nGroup))
                        {
                            // Store each such intersection of complement contributions
                            intersection.push_back(ita->second.vCohoms[i]);
                            intersection.push_back(ita2->second.vCohoms[j]);
                        }
                    }
                }

                // Check if the intersection was able to uniquely identify a contribution cohomology group
                if (intersection.size() == 2)
                {
                    // Apparently, we could uniquely resolve the issue

                    // Add the "normal" monomial to the list of unique ones
                    UniqueContribData newunique;
                    newunique.Cohom.nGroup = intersection[0].nGroup;
                    newunique.Cohom.nFactor = intersection[0].nFactor;
                    newunique.nRationals = ita->second.nRationals;
                    newunique.nRationalsDual = ita->second.nRationalsDual;
                    pair<map<uint64_t, UniqueContribData>::iterator, bool> ret = monoms.unique_monoms.insert(std::pair<uint64_t, UniqueContribData>(monomial, newunique));

                    // Add the "complement" partner monomial to the list of unique ones
                    newunique.Cohom.nGroup = intersection[1].nGroup;
                    newunique.Cohom.nFactor = intersection[1].nFactor;
                    newunique.nRationals = ita2->second.nRationals;
                    newunique.nRationalsDual = ita2->second.nRationalsDual;
                    pair<map<uint64_t, UniqueContribData>::iterator, bool> ret2 = monoms.unique_monoms.insert(pair<uint64_t, UniqueContribData>(complement, newunique));

                    if ((ret.second == false) || (ret2.second == false))
                    {
                        ERR_OUT("INTERNAL: Could not add resolved ambiguous monomials to unique monomials list.");
                        return false;
                    }

                    // Delete both ambiguous cases from the ambiguous list & start the scan from the beginning
                    monoms.ambiguous_monoms.erase(monomial);
                    monoms.ambiguous_monoms.erase(complement);
                    ita = monoms.ambiguous_monoms.begin();
                }
                else
                    ita++;
            }
			else
			{
				// CONSOLE_MSG_OUT("WARNING: The Serre duality check could not entirely resolve all monomials.");
				ita++;
			}
        }
    }

    // In the best case scenario, we could resolve all ambiguous cases and turn them into
    // uniques. Unfortunately, this is not always the case.
    size_t num_ambiguous_after = monoms.ambiguous_monoms.size();
    if (num_ambiguous_after > 0)
    {
        WARN_OUT("The Serre dualization reduction was unable to uniquely resolve " << num_ambiguous_after << " of the original " << num_ambiguous_before << " ambiguous monoms.");
    }

    return true;
}
