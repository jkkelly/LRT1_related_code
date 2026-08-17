#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,asin,sin,log
import numpy as np
from scipy.stats import multivariate_normal,chi2


out0 = open("LRT1_results.txt","w")

#### main
GroupMembers={}
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
    elif cols[0]=="FactorLevels":
        FactorLevels=int(cols[1])
        for j in range(FactorLevels):
            GroupMembers[j]=[]
    elif cols[0]=="Groups":
        for j in range(1,FactorLevels+1):        
            vv=cols[j].split(",")
            for k in range(len(vv)):
                GroupMembers[j-1].append(int(vv[k]))

in1.close()

nullvar={}
in1=open("Nullvar.txt","r")
for line_idx, line in enumerate(in1):
    cols = line.replace('\n', '').split('\t')
# 2	 0.003430499115489819
# 3	 0.0020653571241520574
    nullvar[int(cols[0])]=float(cols[1])
in1.close()


cmat=[] # fill in digonal below for each SNP
for j in range(NoPops):
        row=[ 0.0 for x in range(NoPops) ]
        cmat.append(row) 

out0.write("chrom\tpos")
out0.write('\tNull_model_z0\tNull_model_LL')
out0.write('\tAlt_z\tAlt_LL\tLRT\tpvalue\n')
            
snpcc = 0
src  =open(infile, "r") #
for line_idx, line in enumerate(src):
    cols = line.replace('\n', '').split('\t')

# chr2L	6353	153,123	351,382	199,210	205,186	277,368	142,169	162,215	359,396
    data=[]
    px=[]
    m={}
    p={}
    z={}

    for j in range(2,len(cols)):
        if j in P_in_seq:
            r,a=cols[j].split(",") # 540,46 
            data.append([int(r),int(a)])
            m[j] = int(r)+int(a)
            if m[j]>0:            
                px.append( float(int(r))/float(m[j]) )
                p[j] = float(int(r))/float(m[j])
                z[j] = 2.0*asin(sqrt(p[j]))
                
                
    if len(px)==NoPops: # all pops must have data
        if min(px)<1-minMAF and max(px)>minMAF: # perform test
            snpcc+=1
            out0.write(cols[0]+"\t"+cols[1])
            # Determine MLE for mu under null (mu_0)
            w={}
            sumw=0.0
            z0=0.0 # MLE for common mean
            for j in m:
                var = nullvar[j]+1.0/float(m[j])
                sumw+=(1.0/var)
            for j in m:
                w[j] = (1.0/(nullvar[j]+1.0/float(m[j]))) /sumw
                z0+= w[j]*z[j]

            zvec=[]
            mvec=[]
            cc=0
            for j in m:
                zvec.append(z[j])
                mvec.append(z0)
                cmat[cc][cc]= (nullvar[j]+1.0/float(m[j]))
                cc+=1                         
                     
            LL0 = log(multivariate_normal.pdf(zvec, mvec, cmat))          
            #print(np.average(zvec),z0,LL0)
            out0.write('\t'+str(z0)+"\t"+str(LL0))

            # Determine MLE for mu in each group (mu_{pop})
            LL1=0.0
            z0={}
            zvec=[]
            mvec=[]
            cc=0
            for grp in range(FactorLevels):
                z0[grp]=0.0 # MLE for common mean
                w={}
                sumw=0.0
                for j in GroupMembers[grp]:
                    var = nullvar[j]+1.0/float(m[j])
                    sumw+=(1.0/var)
                for j in GroupMembers[grp]:
                    w[j] = (1.0/(nullvar[j]+1.0/float(m[j]))) /sumw
                    z0[grp]+= w[j]*z[j]    
                for j in GroupMembers[grp]:
                    zvec.append(z[j])
                    mvec.append(z0[grp])
                    cmat[cc][cc]= (nullvar[j]+1.0/float(m[j]))
                    cc+=1                       
                    
            LL1 = log(multivariate_normal.pdf(zvec, mvec, cmat))          
            # print(z0,LL1)
            LRT=2*(LL1-LL0)
            pvalue = 1.0-chi2.cdf(LRT, FactorLevels-1)
            out0.write('\t'+str(z0)+"\t"+str(LL1)+'\t'+str(LRT)+'\t'+str(pvalue)+'\n')
            
src.close()


print("tested snps",snpcc)
out0.close()

