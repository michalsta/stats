# Copyright (c) 2016, Michal Startek
# 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
#    * Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice,
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

__version__ = '1.0'


import math
from collections import Counter
from random import uniform
from numpy.random import binomial


def _beta_1_b(b):
    '''Returns a random variate from beta(1, b) distribution
    using the inverse CDF method'''
    return 1.0-uniform(0.0, 1.0)**(1.0/b)

def _safe_binom(n, p):
    '''Draw a sample from binomial distribution. Doesn't crash
    if p > 1.0 due to numerical inaccuracies.'''
    if p >= 1.0:
        return n
    return binomial(n, p)

def sample_with_replacement(population, probabilities, sample_size):
    '''Performs sampling with replacement from population argument, with
    associated probabilities from second argument. The probabilities must 
    sum to 1. Returns a Counter defining how many times each element has 
    been picked.'''
    pprob = 0.0
    cprob = 0.0
    pidx = 0
    ret = Counter()
    while sample_size > 0:
        pprob += probabilities[pidx]
        # Beta mode
        while (pprob - cprob) * sample_size / (1.0 - cprob) < 1.0:
            cprob += _beta_1_b(sample_size) * (1.0 - cprob)
            while pprob < cprob:
                pidx += 1
                pprob += probabilities[pidx]
            ret[population[pidx]] += 1
            sample_size -= 1
            if sample_size == 0: break
        if sample_size == 0: break
        # Binomial mode
        nrtaken = _safe_binom(sample_size, (pprob-cprob)/(1.0-cprob))
        if nrtaken > 0:
            ret[population[pidx]] += nrtaken
            sample_size -= nrtaken
        pidx += 1
        cprob = pprob
    return ret
    



def sample_with_replacement_online(population, probabilities, sample_size):
    '''Performs sampling with replacement from population argument, with
    associated probabilities from second argument. The probabilities must 
    sum to 1. Yields a stream of tuples: (population_member, times_chosen).
    Accepts generators as first and second argument.
    '''
    p = None
    number = 0
    for popm, times in _sample_with_replacement_online_impl(population, probabilities, sample_size):
        if p is popm:
            number += times
        else:
            if number > 0:
                yield (p, number)
            p = popm
            number = times
    if number > 0:
        yield (p, number)


def _sample_with_replacement_online_impl(population, probabilities, sample_size):
    '''Performs sampling with replacement from population argument, with
    associated probabilities from second argument. The probabilities must 
    sum to 1. Yields a stream of tuples: (population_member, times_chosen).
    Accepts generators as first and second argument. May return duplicate
    tuples and tuples with times_chosen == 0.
    '''
    pprob = 0.0
    cprob = 0.0
    pidx = 0
    population_iter = population.__iter__()
    probabilities_iter = probabilities.__iter__()
    population_next = next(population_iter)
    while sample_size > 0:
        pprob += next(probabilities_iter)
        # Beta mode
        while (pprob - cprob) * sample_size / (1.0 - cprob) < 1.0:
            cprob += _beta_1_b(sample_size) * (1.0 - cprob)
            while pprob < cprob:
                population_next = next(population_iter)
                pprob += next(probabilities_iter)
            yield (population_next, 1)
            sample_size -= 1
            if sample_size == 0: break
        if sample_size == 0: break
        # Binomial mode
        nrtaken = _safe_binom(sample_size, (pprob-cprob)/(1.0-cprob))
        if nrtaken > 0:
            yield (population_next, nrtaken)
            sample_size -= nrtaken
        population_next = next(population_iter)
        cprob = pprob
    

