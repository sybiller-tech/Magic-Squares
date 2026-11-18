#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 13:16:07 2026

@author: sybille

GOAL: Naive computation of the Erhrart polynomial Hn(r) associated with magic squares
"""
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import comb
import time



import sympy as sp

n =7

x = sp.symbols("x")



def rows_with_sum(r, caps):
    """generates every way to fill one row of the matrix with line sums R, given upper rows of column sums R-CAPS:"""
    """--> yields all tuples of nonnegative integers that sum to r, where entry i can’t exceed caps[i] (the amount still available in column i)."""
    n = len(caps)
    def rec(i, remaining):
        # decides entries i, i+1, ..., n-1 so that they sum to remaining
        if i == n:
            if remaining == 0: #Otherwise yield nothing
                yield ()
            return
        later_cap = sum(caps[i + 1:]) # remaining capacity of columns after column i
        lo = max(0, remaining - later_cap) #v must take up at least remaining - later_cap
        hi = min(caps[i], remaining) #v can’t exceed the column’s cap, and can’t exceed what’s left to distribute.
        for v in range(lo, hi + 1):
            for rest in rec(i + 1, remaining - v):
                yield (v,) + rest
    yield from rec(0, r) #Start at position 0 with the full sum r, and pass every result up to the caller


def H(n, r):
    """Number of n x n nonnegative integer matrices with all row/column sums r."""
    if n == 0:
        return 1
    if r == 0:
        return 1

    @lru_cache(maxsize=None)
    def count(rows_left, cols):
        # cols: sorted tuple of remaining column sums
        if rows_left == 1:
            return 1  # last row is forced (sum(cols) == r by invariant)
        total = 0
        for row in rows_with_sum(r, cols):
            new_cols = tuple(sorted(c - x for c, x in zip(cols, row)))
            total += count(rows_left - 1, new_cols)
        return total

    return count(n, tuple([r] * n))

def h(n):
    """h-vector of the polytope K_n"""
    h=[]
    D= n**2-3*n+2 #length of the h-vector of K_n -1
    N=(n-1)**2+1
    print("D=",D)
    start = time.process_time()
    Hvals = [H(n, r) for r in range(0,D//2+1)]   
    result = sum(range(10_000_000))
    end = time.process_time()
    print(f"Elapsed: {end - start:.4f} seconds")
    for r in range(D//2+1):
        #print("r=", r)
        h.append(sum((-1)**k * comb(N, k) * Hvals[r - k] for k in range(r + 1)))
    for i in range(D//2,D):
        #print(i)
        h.append(h[D-i-1])
    print("Delta=")
    for i in range(1,D-1):
        Delta= h[i]**2-h[i-1]*h[i+1]
        #print(Delta)
        assert(Delta >=0)
    return h



if __name__ == "__main__":
    # Table of values
    for n in range(1, n+1):
        hvect=h(n)
        print(n)
        #Hpol=sp.interpolate(Hvals, x)
        #print(f"H_{n}(r), r=0..{d+1}:", Hvals)
        print(f"h_{n}:", hvect)
        print(f"Volume of K_{n}:", sum(hvect[i] for i in range(n**2-3*n+2+1)))

        

    # # Sanity check against the known closed forms for n =2, 3:
    assert all(H(2, r) == r+1 for r in range(10))
    assert all(H(3, r) == comb(r + 2, 4) + comb(r + 3, 4) + comb(r+4,4) for r in range(10))



    
