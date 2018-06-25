import numpy as np
import gmres


def f_traction(A,B, ux, uy):
    n = len(ux)
    u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    u[:n,0], u[n:,0] = ux[:], uy[:]
    
    b = np.matmul(B,u)
    f = np.linalg.solve(A, b)
    
    #print b
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n,0], f[n:,0]
    
    #print fx
    #print
    #print fy
    
    return fx, fy
    
def f_traction_gmres(A,B, ux, uy):
    n = len(ux)
    nn = 2*n
    u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    u[:n,0], u[n:,0] = ux[:], uy[:]
    f = np.zeros(nn, dtype = np.float64)
    f[:] = 1.0E-16
    
    b = np.matmul(B,u)
   
    
    #D = np.zeros([nn,nn], dtype=np.float64)
    #for i in range(nn):
    #    D[i,i] = A[i,i]**(-1)


    #AA = np.matmul(D,A)
    #bb = np.matmul(D,b)
    b2 = np.transpose(b)
    
    #print b
    #print b2
    f,_ = gmres_mgs(A,f,b2, nn, 500, 1.0E-14)
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n], f[n:]
    
    del f
    
    #print fx
    #print
    #print fy
    
    return fx, fy
    
def tangential_of_ft(fx,fy, tx, ty, h):
    ftn = np.zeros(len(fx), dtype = np.float64)
    dDrag = np.zeros(len(fx), dtype = np.float64)    
    
    ftn[:] = fx[:]*tx[:] + fy[:]*ty[:]
#    print fx
#    print 
#    print fy
#    print 
#    print ftn    
    
    dDrag[:] = fx[:]*h[:]
    
    return sum(dDrag)
    
def calc_ndrag(fx,h):
    dDrag = np.zeros(len(fx), dtype = np.float64)
    dDrag[:] = fx[:]*h[:]
    
    return sum(dDrag)
    