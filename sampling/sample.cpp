 /* Copyright (c) 2016, Michal Startek
  * 
  * All rights reserved.
  * 
  * Redistribution and use in source and binary forms, with or without modification,
  * are permitted provided that the following conditions are met:
  * 
  *    * Redistributions of source code must retain the above copyright notice,
  *      this list of conditions and the following disclaimer.
  *    * Redistributions in binary form must reproduce the above copyright notice,
  *      this list of conditions and the following disclaimer in the documentation
  *      and/or other materials provided with the distribution.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
  * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  */



#include <random>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <unordered_map>

std::random_device rd;
std::mt19937 gen(rd());
boost::mt19937 rng;

std::uniform_real_distribution<double> stdunif(0.0, 1.0);


inline double beta_1_b(double b)
{
	return 1.0 - pow(stdunif(gen), 1.0/b);
}


inline int binom(int tries, double succ_prob)
{
	if (succ_prob >= 1.0)
		return tries;
	boost::random::binomial_distribution<> bd(tries, succ_prob);
	boost::variate_generator<boost::mt19937&, boost::random::binomial_distribution<> > var_b(rng, bd);
	return var_b();
}





void sample(const double* probs, int* sample, int samples, double switchover = 1.0)
/* Generates a random sample from population 0..n-1 with probabilities probs, and saves it to
 * table sample (which must be of size equal to samples parameter). Parameter switchover controls 
 * preference between beta and binomial modes.
 */
{
	double pprob = 0.0;
	double cprob = 0.0;
	unsigned int pidx = 0;
	unsigned int sampleidx = 0;
	while(samples > 0)
	{
		pprob += probs[pidx];
		while(((pprob - cprob) * samples / (1.0 - cprob)) < switchover)
		{
			cprob += beta_1_b(samples) * (1.0 - cprob);
			while(pprob < cprob)
				pprob += probs[++pidx];
			sample[sampleidx++] = pidx;
			samples--;
			if(samples == 0)
				break;
		}
		if(samples == 0)
			break;
		unsigned int nrtaken = binom(samples, (pprob-cprob)/(1.0-cprob));
		for(unsigned int i=0; i<nrtaken; i++)
			sample[sampleidx++] = pidx;
		samples -= nrtaken;
		pidx++;
		cprob = pprob;
	}
}

void sample_cntr(const double* probs, int* sample, int popsize, int samples, double switchover)
/* Generates a random sample from population 0..n-1 with probabilities probs, and saves counts to
 * the table sample (which must be of size popsize). Parameter switchover controls
 * preference between beta and binomial modes.
 */
{
	memset(sample, 0, sizeof(int)*popsize);
        double pprob = 0.0;
        double cprob = 0.0;
        unsigned int pidx = 0;
        while(samples > 0)
        {
                pprob += probs[pidx];
                while(((pprob - cprob) * samples / (1.0 - cprob)) < switchover)
                {
                        cprob += beta_1_b(samples) * (1.0 - cprob);
                        while(pprob < cprob)
                                pprob += probs[++pidx];
			sample[pidx] += 1;
                        samples--;
                        if(samples == 0)
                                break;
                }
                if(samples == 0)
                        break;
                unsigned int nrtaken = binom(samples, (pprob-cprob)/(1.0-cprob));
		sample[pidx] += nrtaken;
                samples -= nrtaken;
                pidx++;
                cprob = pprob;
        }
}

std::unordered_map<int, int>* sample_ht(const double* probs, int samples, double switchover)
/* Generates a random sample from population 0..n-1 with probabilities probs, and returns a hashtable
 * containing the counts. Parameter switchover controls preference between beta and binomial modes.
 */
{
	std::unordered_map<int, int>* sampl_multiset = new std::unordered_map<int, int>();
        double pprob = 0.0;
        double cprob = 0.0;
        unsigned int pidx = 0;
        while(samples > 0)
        {
                pprob += probs[pidx];
                while(((pprob - cprob) * samples / (1.0 - cprob)) < switchover)
                {
                        cprob += beta_1_b(samples) * (1.0 - cprob);
                        while(pprob < cprob)
                                pprob += probs[++pidx];
			if(sampl_multiset->find(pidx) == sampl_multiset->end())
				(*sampl_multiset)[pidx] = 1;
			else
				(*sampl_multiset)[pidx] += 1;
                        samples--;
                        if(samples == 0)
                                break;
                }
                if(samples == 0)
                        break;
                unsigned int nrtaken = binom(samples, (pprob-cprob)/(1.0-cprob));
		if(nrtaken > 0)
		{
			if(sampl_multiset->find(pidx) == sampl_multiset->end())
				(*sampl_multiset)[pidx] = nrtaken;
			else
				(*sampl_multiset)[pidx] += nrtaken;
		}
                samples -= nrtaken;
                pidx++;
                cprob = pprob;
        }
        return sampl_multiset;

}

