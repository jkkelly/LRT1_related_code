#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from math import sqrt,asin,sin
import numpy as np




def perform_ols(incidence_matrix: np.ndarray, y_values: np.ndarray):
    """Performs Ordinary Least Squares (OLS) regression using the normal equation:

    beta = (X^T * X)^(-1) * X^T * y

    Parameters:
     incidence_matrix (np.ndarray): Design/Incidence matrix X of shape (N, K)
     y_values (np.ndarray): Observed vector y of shape (N,) or (N, 1)

    Returns:
     dict: Dictionary containing estimated coefficients, fitted values, residuals, and RSS.
    """
    X = np.array(incidence_matrix, dtype=float)
    y = np.array(y_values, dtype=float).reshape(-1, 1)

    # Validate dimensions
    if X.shape[0] != y.shape[0]:
        raise ValueError(
            f"Dimension mismatch: X has {X.shape[0]} rows, but y has {y.shape[0]} values."
        )

    # Check for full column rank to ensure invertibility
    rank = np.linalg.matrix_rank(X)
    if rank < X.shape[1]:
        raise np.linalg.LinAlgError(
            f"Matrix X is rank deficient (Rank = {rank}, Columns = {X.shape[1]}). "
            "Cannot invert (X^T * X)."
        )

    # Compute parameters using least squares solver (more numerically stable than direct inverse)
    beta, rss, rank, s = np.linalg.lstsq(X, y, rcond=None)

    # Compute fitted values and residuals
    y_hat = X @ beta
    residuals = y - y_hat

    return {
        "beta": beta.flatten(),
        "fitted_values": y_hat.flatten(),
        "residuals": residuals.flatten(),
        "rss": float(np.sum(residuals**2)),
    }


#### main

in1=open("Control.txt")
for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
    if cols[0]=="infile":
        infile = cols[1] # "example.data.txt" 
    elif cols[0]=="minMAF":
        minMAF=float(cols[1])
    elif cols[0]=="NoPops":
        NoPops=int(cols[1])
    elif cols[0]=="PCols":
        P_in_seq=[]
        vv=cols[1].split(",")
        if len(vv) != NoPops:
            print("error : populations")
            break
        for j in range(len(vv)):
            P_in_seq.append(int(vv[j]))

    # print(cols[0],"heynow")
in1.close()
out2=open("Nullvar.txt","w")
# out1=open("pairdz.byfreq.txt","w")



dz={}
dz2={} # squared divergence
rdv={} # read depth var
for j in range(NoPops-1):
    for k in range(j+1,NoPops):
        dz2[str(j)+"_"+str(k)]=0.0
        rdv[str(j)+"_"+str(k)]=0.0        
        dz[str(j)+"_"+str(k)] =[]

snpcc = 0
src  =open(infile, "r") #
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

# chr2L	6353	153,123	351,382	199,210	205,186	277,368	142,169	162,215	359,396
    data=[]
    px=[]
    for j in range(2,len(cols)):
        if j in P_in_seq:
            r,a=cols[j].split(",") # 540,46 
            data.append([int(r),int(a)])
            m0 = int(r)+int(a)
            if m0>0:            
                px.append( float(int(r))/float(m0) )
            
    if len(px)==NoPops: # all pops must have data
        pm=np.average(px)
        if min(pm,1-pm)>minMAF:
            snpcc+=1
            for j in range(NoPops-1):
                m0 = sum(data[j])
                p0 = float(data[j][0])/float(m0)
                z0 = 2.0*asin(sqrt(p0))
                for k in range(j+1,NoPops):
                    m1=sum(data[k])
                    p1 = float(data[k][0])/float(m1)
                    z1 = 2.0*asin(sqrt(p1))
                    if line_idx % 2 == 0:
                        dz[str(j)+"_"+str(k)].append(z1-z0)
                    else:
                        dz[str(j)+"_"+str(k)].append(z0-z1)                        
                    dz2[str(j)+"_"+str(k)]+= (z1-z0)**2.0
                    rdv[str(j)+"_"+str(k)]+= 1.0/float(m1)+1.0/float(m0)
            
src.close()

yx=[]
xx=[]
print("included snps",snpcc)
for j in range(NoPops-1):
    for k in range(j+1,NoPops):
        xx.append([0 for x in range(NoPops)])
        xx[-1][j]=1
        xx[-1][k]=1        
        dz2[str(j)+"_"+str(k)]= dz2[str(j)+"_"+str(k)]/float(snpcc)
        rdv[str(j)+"_"+str(k)]= rdv[str(j)+"_"+str(k)]/float(snpcc)                    

        ranked_z=sorted(dz[str(j)+"_"+str(k)])
        n25=ranked_z[int(snpcc/4)]
        n50=ranked_z[int(snpcc/2)]
        n75=ranked_z[int(3*snpcc/4)]
        yx.append( ((n75-n25)/1.349)**2.0 - rdv[str(j)+"_"+str(k)] ) # esimate of v(i) + v(j)
                    
        print (j,k,"mean dz2",dz2[str(j)+"_"+str(k)],"var_IQR ",((n75-n25)/1.349)**2.0,"readdepth_var",rdv[str(j)+"_"+str(k)])

        
# Run OLS calculation
try:
    results = perform_ols(np.array(xx), np.array(yx))

    print("\n--- Results ---")
    print("Estimated null vars:")

    for j in range(NoPops):
        out2.write(str(P_in_seq[j])+'\t'+str(results["beta"][j])+'\n')

    for i, b in enumerate(results["beta"]):
        print(f"  beta_{i} = {b:.6f}")



except Exception as e:
    print(f"\nError computing OLS: {e}")
        
        
        
out2.close()

