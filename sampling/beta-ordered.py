from random import uniform
from math import exp, log

try:
    xrange
except NameError:
    # We're in Python3
    xrange = range


def ordered_unif_seq(samplesize):
    if samplesize <= 0:
        return
    last = uniform(0.0, 1.0)**(1.0/samplesize)
    yield 1.0-last
    for n in xrange(samplesize-1, 0, -1):
        last = last * uniform(0.0, 1.0)**(1.0/n)
        yield 1.0-last


def ordered_unif_exp(samplesize):
    if samplesize <= 0:
        return
    L = [-log(uniform(0.0, 1.0))]
    for _unused in xrange(samplesize):
        L.append(L[-1] - log(uniform(0.0, 1.0)))
    last = L.pop()
    return [x/last for x in L]



def _beta_1_b(b):
    '''Returns a random variate from beta(1, b) distribution
    using the inverse CDF method'''
    return 1.0-uniform(0.0, 1.0)**(1.0/b)


def ordered_unif_beta(samplesize):
    now = 0.0
    for n in xrange(samplesize-1, 0, -1):
        B = uniform(0.0, 1.0)**(1.0/n)
        now = 1.0 - B*(1.0 - now)
        yield now

def ordered_unif_seq_rev(samplesize):
    if samplesize <= 0:
        return
    last = uniform(0.0, 1.0)**(1.0/samplesize)
    yield last
    for n in xrange(samplesize-1, 0, -1):
        last = last * uniform(0.0, 1.0)**(1.0/n)
        yield last


print(list(ordered_unif_seq(100)))
