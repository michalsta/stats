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
    yield last
    for n in xrange(samplesize-1, 0, -1):
        last = last * uniform(0.0, 1.0)**(1.0/n)
        yield last


def ordered_unif_exp(samplesize):
    if samplesize <= 0:
        return
    csum = -log(uniform(0.0, 1.0))
    for _unused in xrange(samplesize):
        expv = -log(uniform(0.0, 1.0))
        next_csum = csum + expv
        print(csum, next_csum, csum/next_csum)
        yield csum/next_csum
        csum = next_csum

print(list(ordered_unif_exp(100)))
