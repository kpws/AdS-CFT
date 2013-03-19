import numpy as np
cimport numpy as np


cdef extern from "solveODE.h":
    int solveODE(double * start, double eps, int n, double *out, double *params)

cdef extern from "stdlib.h":
        ctypedef int size_t
        void* malloc(size_t)
        void free(void*)

def solve(start, eps, n, params):
    dim=len(start)
    cdef double* out_v = <double*> malloc(dim*n*sizeof(double))
    cdef double* start_v = <double*> malloc(dim*sizeof(double))
    cdef double* params_v = <double*> malloc(len(params)*sizeof(double))
    cdef int i=0
    for ii in start:
        start_v[i]=<double> ii
        i+=1
    i=0
    for ii in params:
        params_v[i]=<double> ii
        i+=1

    solveODE(start_v, eps, n, out_v, params_v)
    ret=np.linspace(1-eps,eps,n), [[out_v[i*dim+j] for j in range(dim)] for i in range(n)]
    free(out_v)
    free(start_v)
    return ret
