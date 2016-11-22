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
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <unordered_map>


#define DO_CLOCK(str) { std::cout << str << ": " << ((double)(clock() - timer))/CLOCKS_PER_SEC << std::endl; timer = clock(); }

std::random_device rd;
std::mt19937 gen(rd());
boost::random::mt19937 rng(rd());

std::uniform_real_distribution<double> stdunif(0.0, 1.0);


inline double beta_1_b(double b, std::mt19937 rdev = gen)
{
	return 1.0 - pow(stdunif(rdev), 1.0/b);
}


inline int binom(int tries, double succ_prob, boost::random::mt19937 rdev = rng)
{
	if (succ_prob >= 1.0)
		return tries;
	boost::random::binomial_distribution<> bd(tries, succ_prob);
	boost::variate_generator<boost::mt19937&, boost::random::binomial_distribution<> > var_b(rdev, bd);
	return var_b();
}





void sample(const double* probs, unsigned int* samplespace, unsigned int samples, double switchover = 1.0)
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
			samplespace[sampleidx++] = pidx;
			samples--;
			if(samples == 0)
				break;
		}
		if(samples == 0)
			break;
		unsigned int nrtaken = binom(samples, (pprob-cprob)/(1.0-cprob));
		for(unsigned int i=0; i<nrtaken; i++)
			samplespace[sampleidx++] = pidx;
		samples -= nrtaken;
		pidx++;
		cprob = pprob;
	}
}

void sample_cntr(const double* probs, unsigned int* sample, unsigned int popsize, unsigned int samples, double switchover = 1.0)
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

std::unordered_map<int, int>* sample_ht(const double* probs, unsigned int samples, double switchover)
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


/* 
 * The following is a multithreaded implementation of the above
 */




#ifdef THREADS
#include <pthread.h>


struct threadargs_t
{
        const double* probs;
	unsigned int probsstart;
	unsigned int probsend;
        unsigned int* samplespace;
        unsigned int totalsamples;
        double switchover;
        double* fragment_probs;
	unsigned int* fragment_sizes;
        unsigned int thread_id;
	pthread_barrier_t* barrier;
	unsigned int no_threads;
	unsigned int pidx;
};

template<typename T> void print_array(T* array, unsigned int size)
{
	std::cout << "-------------------------------------" << std::endl;
	for(int ii=0; ii<size; ii++)
		std::cout << array[ii] << " ";
	std::cout << std::endl << "-------------------------------------" << std::endl;
}


void sample_threadfunc(threadargs_t* args)
{
	clock_t timer = clock();
	std::random_device rd;
	std::mt19937 gen(rd());
	boost::random::mt19937 rng(rd());

	double myprobs = 0.0;
	for (int ii = args->probsstart; ii < args->probsend; ii++)
		myprobs += args->probs[ii];

	if(args->thread_id == 0)
		DO_CLOCK("After MT init")
	
	args->fragment_probs[args->thread_id] = myprobs;

	int leader = pthread_barrier_wait(args->barrier);

	if(args->thread_id == 0)
		DO_CLOCK("First barrier");

	if(leader == PTHREAD_BARRIER_SERIAL_THREAD)
		/* Deliberately not using switchover arg. Even if input has weird distribution which justifies it,
		 * this distribution is different. */
		sample_cntr(args->fragment_probs, args->fragment_sizes, args->no_threads, args->totalsamples);

	pthread_barrier_wait(args->barrier);

	if(args->thread_id == 0)
		DO_CLOCK("Sample_cntr, second barrier");

	/* Probably faster to compute this in each thread from scratch, like this,
	 * than to synchronise again just to compute this... */
	unsigned int* samplespace = args->samplespace;
	for(int ii=0; ii<args->thread_id; ii++)
		samplespace += args->fragment_sizes[ii];

	unsigned int samples = args->fragment_sizes[args->thread_id];
	double switchover = args->switchover;
	
	/* Sadly, not *quite* the same as the single-threaded version above... */
	double correction = 1.0/myprobs;
        double pprob = 0.0;
        double cprob = 0.0;
	const double* probs = args->probs;
	
        unsigned int pidx = args->probsstart;
        unsigned int sampleidx = 0;
        while(samples > 0)
        {
                pprob += probs[pidx]*correction;
                while(((pprob - cprob) * samples / (1.0 - cprob)) < switchover)
                {
                        cprob += beta_1_b(samples, gen) * (1.0 - cprob);
                        while(pprob < cprob)
                                pprob += probs[++pidx]*correction;
                        samplespace[sampleidx++] = pidx;
                        samples--;
                        if(samples == 0)
                                break;
                }
                if(samples == 0)
                        break;
                unsigned int nrtaken = binom(samples, (pprob-cprob)/(1.0-cprob), rng);
                for(unsigned int i=0; i<nrtaken; i++)
                        samplespace[sampleidx++] = pidx;
                samples -= nrtaken;
                pidx++;
                cprob = pprob;
        }

	if(args->thread_id == 0)
	        DO_CLOCK("After main loop");

}

void sample_multithreaded(const double* probs, unsigned int* samplespace, unsigned int popsize, unsigned int samples, double switchover = 1.0, unsigned int no_threads = 1)
{
	
        if(no_threads == 1)
                return sample(probs, samplespace, samples, switchover);

        threadargs_t* threadargs = new threadargs_t[no_threads];
        double* fragment_probs = new double[no_threads];
	unsigned int* fragment_sizes = new unsigned int[no_threads];

        pthread_barrier_t barrier;
        pthread_barrier_init(&barrier, nullptr, no_threads);

        pthread_t* thread_ids = new pthread_t[no_threads];

	unsigned int perthread = popsize / no_threads;
	if(popsize % no_threads > 0)
		perthread++;

        for(unsigned int ii=0; ii<no_threads; ii++)
        {
                threadargs[ii].thread_id = ii;
		threadargs[ii].probs = probs;
		threadargs[ii].probsstart = perthread*ii;
		threadargs[ii].probsend = std::min(perthread*(ii+1), popsize);
		threadargs[ii].samplespace = samplespace;
		threadargs[ii].totalsamples = samples;
		threadargs[ii].switchover = 1.0;
		threadargs[ii].fragment_probs = fragment_probs;
		threadargs[ii].fragment_sizes = fragment_sizes;
		threadargs[ii].barrier = &barrier;
		threadargs[ii].no_threads = no_threads;

		assert(pthread_create(&thread_ids[ii], NULL, (void* (*)(void*)) sample_threadfunc, &threadargs[ii]) == 0);

        }

	for(unsigned int ii=0; ii<no_threads; ii++)
	{
		pthread_join(thread_ids[ii], nullptr);
	}


	delete[] threadargs;
	delete[] fragment_probs;
	delete[] fragment_sizes;
	delete[] thread_ids;
	pthread_barrier_destroy(&barrier);

}
#endif /* THREADS */




int main(int argc, char** argv)
{
	int popsize = atoi(argv[1]);
	int samplesize = atoi(argv[2]);
	int threads = atoi(argv[3]);

	double* probs = new double[popsize];
	unsigned int* samplespace = new unsigned int[samplesize];

	double csum = 0.0;
	for (int ii=0; ii<popsize; ii++)
	{
		probs[ii] = stdunif(gen);
		csum += probs[ii];
	}
	for (int ii=0; ii<popsize; ii++)
		probs[ii] /= csum;

	clock_t timer = clock();
	sample_multithreaded(probs, samplespace, popsize, samplesize, 1.0, threads);
	std::cout << ((double)(clock() - timer))/CLOCKS_PER_SEC << std::endl;

//	for (int ii=0; ii<samplesize; ii++)
//		std::cout << samplespace[ii] << std::endl;

	delete[] probs;
	delete[] samplespace;
	
}
