import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

na=50
neta=2
neps=2
neg=-1e10
sigma=2
gam=0.31
t=45
tr=25
nage=t+tr
taur=0.429
assetmax=10
assetmin=0

agrid=np.linspace(0,assetmax,na)

def u(x,h):
    return ((x)**gam*(1-h)**(1-gam))**(1-sigma)/(1-sigma)
		

def getoptpolicy(wseq,rseq,penseq,trseq,tauwseq,taubseq):
    
    copt0=np.zeros((nage,neta,neps,na))
    val0=np.zeros((nage,neta,neps,na))
    
    for ia in range(na):
        for ieta in range(neta):
            for ieps in range(neps):
                c=(1+(1-taur)*rseq[nage])*agrid[ia]+penseq[nage,ieps]+trseq[nage]
                copt0[nage,ieta,ieps,ia]=c	
                y=u(c,0)
                if c>0:
                    val0[nage,ieta,ieps,ia]=y
                else:
                    val0[nage,ieta,ieps,ia]=neg
    
    for it in range(tr,0,-1):
        for ia in range(na):
            for ieta in range(neta):
                for ieps in range(neps):
                    cmax=(1+(1-taur)*rseq[it+t])*agrid[ia]+penseq[it+t,ieps]+trseq[it+t]-assetmin
                    
                    if cmax>0:
                        m0=0
                    vr11=val0[t+it,ieta,ieps,:]
                    df=vr11[1:na-1]-vr11[0:na-2]
                    
                    if sum(df<0)>0:
                        print('Not Monotone')
                    
        
    
