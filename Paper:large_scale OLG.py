import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate 
import os


data=pd.read_excel('exp.xlsx')
sp1=pd.read_excel('sp.xlsx')['sp']
efage=pd.read_excel('efage.xlsx')['efage']

psi=0.001
alpha=0.36
delta=0.08
na=50
nin=1000
neta=2
neps=2
neg=-1e6
sigma=2
gam=0.31
t=45
tr=25
nage=t+tr
taur=0.429
assetmax=10
pen_repl=0.5
gn=0.02
#na1 is the zero position
na1=0

garate=0.02
beta1=1.011
pi_eta=np.array([[0.98,0.02],[0.02,0.98]])
eps=np.array([0.57,1.43])
eta=np.array([0.727,1.273])



wseq=data['w']
rseq=data['r']
penseq=np.array((data['pen1'],data['pen2'])).T
trseq=data['tr']
tauwseq=data['tauw']
taubseq=data['taub']
agrid=np.linspace(0.0001,assetmax,na)
incomemax=4
incomegrid=np.linspace(0,incomemax,nin)
assetmin=agrid[0]

mass=np.ones((nage,1))
for i in range(1,nage):
    mass[i]=mass[i-1]*sp1[i-1]/(1+gn)
mass=mass/(sum(mass))
massinitial=mass

def golden(f,ay,by,cy,tol=1e-5):
    r1=2/(5**0.5+1)
    r2=1-r1
    x0=ay
    x3=cy
    if abs(cy-by)<=abs(by-ay):
        x1=by
        x2=by+r2*(cy-by)
    else:
        x2=by
        x1=by-r2*(by-ay)
    f1=-f(x1)
    f2=-f(x2)
    
    while abs(x3-x0)>tol*max(1,abs(x1)+abs(x2)):
        if f2<f1:
            x0=x1
            x1=x2
            x2=r1*x1+r2*x3
            f1=f2
            f2=-f(x2)
        else:
            x3=x2
            x2=x1
            x1=r1*x2+r2*x0
            f2=f1
            f1=-f(x1)
        
    if f1<f2:
        xmin=x1
    else:
        xmin=x2
    return xmin
    
def u(x,h):
    value_u=(x)**gam*(1-h)**(1-gam)
    if value_u<0:
        value_u=0
    return (value_u)**(1-sigma)/(1-sigma)
 
def rvalue(a):
    if a==assetmin:
        return vr11[1]
    if a>=assetmax:
        return vr11[na-1]
    ia0=sum(a>agrid)
    lam=(agrid[ia0]-a) / (agrid[ia0]-agrid[ia0-1])
    if ia0==na:
        expv=vr11[ia0-1]
    else:
        expv=lam*vr11[ia0-1]+(1-lam)*vr11[ia0]
    return expv

def value1(x):
    c=(1+(1-taur)*rseq[it+t-1])*agrid[ia]+penseq[it+t-1,ieps]+trseq[it+t-1]-x*(1+garate)	
    if c<0:
        return neg
    y=u(c,0)+beta1*sp1[t+it-1]*(1+garate)**(gam*(1-sigma))*rvalue(x)
    return y	

def wvalue(a1):
    c=gam*( (1+(1-taur)*rseq[it])*asset0+trseq[it]+(1-tauwseq[it]-taubseq[it])*w0-a1*(1+garate) )
    l0=1-c/ ( (1-tauwseq[it]-taubseq[it])*w0 ) * (1-gam)/gam
    
    if l0<0:
        l0=0
        c= (1+(1-taur)*rseq[it])*asset0+trseq[it]-a1*(1+garate)
    
    #f1=interpolate.interp1d(agrid,v11e,'linear')
    #f2=interpolate.interp1d(agrid,v12e,'linear')
    #x1=f1(a1)
    #x2=f2(a1)
    #'''
    i_a=0
    while i_a<na:
        if a1>agrid[i_a]:
            i_a=i_a+1
        else:
            break
    x1=(a1-agrid[i_a])/(agrid[i_a-1]-agrid[i_a])*v11e[i_a-1]+\
        (a1-agrid[i_a-1])/(agrid[i_a]-agrid[i_a-1])*v11e[i_a]
    x2=(a1-agrid[i_a])/(agrid[i_a-1]-agrid[i_a])*v12e[i_a-1]+\
        (a1-agrid[i_a-1])/(agrid[i_a]-agrid[i_a-1])*v12e[i_a]
    #'''
    expected_value=pi_eta[ieta,0]*x1+pi_eta[ieta,1]*x2
    y=u(c,l0)
    y=y+beta1*sp1[it]*(1+garate)**(gam*(1-sigma))*expected_value
    if c<0:
        y=neg
    return y

def getoptpolicy(wseq,rseq,penseq,trseq,tauwseq,taubseq):
    
    global vr11, v11e, v12e, m, w0, WV, asset0
    global ia,ieta,ieps,it
    
    copt0=np.zeros((nage,neta,neps,na))
    lopt0=np.zeros((t,neta,neps,na))
    val0=np.zeros((nage,neta,neps,na))
    aopt0=np.zeros((nage,neta,neps,na))
    
    for ia in range(na):
        for ieta in range(neta):
            for ieps in range(neps):
                c=(1+(1-taur)*rseq[nage-1])*agrid[ia]+penseq[nage-1,ieps]+trseq[nage-1]
                copt0[nage-1,ieta,ieps,ia]=c	
                y=u(c,0)
                if c>0:
                    val0[nage-1,ieta,ieps,ia]=y
                else:
                    val0[nage-1,ieta,ieps,ia]=neg
    #------------------------------------------------------------------------
                    
    #Part 1 : Retirement Problem
    
    #------------------------------------------------------------------------
    
    for it in range(tr-1,0,-1):
        #print("vr11",vr11)
        print("old",it)
        for ia in range(na):
            for ieta in range(neta):
                for ieps in range(neps):
                    #cmax=(1+(1-taur)*rseq[it+t])*agrid[ia]+penseq[it+t,ieps]+trseq[it+t]-assetmin
                    #if cmax>0:
                    #vr=val0.reshape(((t+tr)*neps*neta,na))  
                    #vr11=vr[(it+t)*neta*neps+(ieta-1)*neps+ieps,:]
                    vr11=val0[t+it,ieta,ieps,:]
                    #print("vr",vr[-10:-1,:])
                    
                    df=vr11[1:na-1]-vr11[0:na-2]
                    
                    if sum(df<0)>0:
                        print('Not Monotone')
                        break
                        
                    ax=0
                    bx=-1
                    cx=-2
                    v0=neg
                    m=0
                    while ax>bx or bx>cx:
                        
                        m=m+1
                        v1=value1(agrid[m])
                        #print(v1)
                        #AX.append(v1)
                        if v1>v0:
                            if m==0:
                                ax=agrid[m]
                                bx=agrid[m]
                            else:
                                bx=agrid[m]
                                ax=agrid[m-1]
                            v0=v1
                            #m0=m
                        else:
                            cx=agrid[m]
                        if m==(na-1):
                            ax=agrid[m-1]
                            bx=agrid[m]
                            cx=agrid[m]
                            break
                    #Judge the position of the m
                    if ax==bx:
                        v1=value1(ax+psi)
                        if v1>v0:
                            bx=ax+psi
                            a1=golden(value1,ax,bx,cx)
                        else:
                            a1=ax
                    elif bx==cx:
                        v0=value1(agrid[na-1])
                        v1=value1(agrid[na-1]-psi)
                        if v0<v1:
                            bx=agrid[na-1]-psi
                            a1=golden(value1,ax,bx,cx)
                        else:
                            a1=agrid[na-1]
                    else:
                        a1=golden(value1,ax,bx,cx)
                        
                    c=(1+(1-taur)*rseq[it+t-1])*agrid[ia]+penseq[it+t-1,ieps]+trseq[it+t-1]-a1*(1+garate)
                    aopt0[it+t-1,ieta,ieps,ia]=a1
                    copt0[it+t-1,ieta,ieps,ia]=c
                    val0[it+t-1,ieta,ieps,ia]=value1(a1)
    #os.system('pause')
    #'''
    #------------------------------------------------------------------------------------
    
    #Part 2: Young's Problem
    
    #------------------------------------------------------------------------------------
    
    for it in range(t-1,-1,-1):
        print("young",it)
        for ia in range(na):
            asset0=agrid[ia]
            for ieps in range(neps):
                v11e=val0[it+1,0,ieps,:]
                v12e=val0[it+1,1,ieps,:]
                
                
                df=v11e[1:na-1]-v11e[0:na-2]
                if sum(df<0)>0:
                    print('v11e Not Monoply')
                df=v12e[1:na-1]-v12e[0:na-2]
                if sum(df<0)>0:
                    print('v12e Not Monoply')
                
                for ieta in range(neta):
                    #print(it,ia,ieps,ieta)
                    w0=wseq[it]*efage[it]*eps[ieps]*eta[ieta]
                    
                    #'''
                    ax=0
                    bx=-1
                    cx=-2
                    vw0=neg
                    m=-1
                    while ax>bx or bx>cx:
                           
                        m=m+1
                        vw1=wvalue(agrid[m])
                        if vw1>vw0:
                            if m==0:
                                ax=agrid[m]
                                bx=agrid[m]
                                
                            else:
                                bx=agrid[m]
                                ax=agrid[m-1]
                            vw0=vw1
                        else:
                            cx=agrid[m]
                        if m==(na-1):
                            ax=agrid[m-1]
                            bx=agrid[m]
                            cx=agrid[m] 
                            break
                    
                    #Check
                    if ax==bx:
                        a1=agrid[0]
                        vw0=wvalue(a1)
                        vw1=wvalue(a1+psi)
                        if vw1>vw0:	
                            bx=assetmin+psi
                            a1=golden(wvalue,ax,bx,cx)
                    elif bx==cx:
                        a1=agrid[na-1]
                        vw0=wvalue(a1)
                        vw1=wvalue(a1-psi)
                        if vw0<vw1 or vw0==vw1:		# corner solution: a'=assetmax?
                            if vw1>neg:		# both solutions for bx and cx imply c<0
                                bx=a1-psi
                            else:
                                bx=ax+psi
                            a1=golden(wvalue,ax,bx,cx)
                    else:
                        a1=golden(wvalue,ax,bx,cx)
                        
                    aopt0[it,ieta,ieps,ia]=a1
                    w0=wseq[it]*efage[it]*eps[ieps]*eta[ieta]

                    c=gam*( (1+(1-taur)*rseq[it])*asset0+trseq[it]+(1-tauwseq[it]-taubseq[it])*w0-a1*(1+garate) )
                    
                    l0=1-c/ ( (1-tauwseq[it]-taubseq[it])*w0 ) * (1-gam)/gam
                    if l0<0:
                        l0=0
                        c= (1+(1-taur)*rseq[it])*asset0+trseq[it]-a1*(1+garate)
                    
                    
                    lopt0[it,ieta,ieps,ia]=l0
                    copt0[it,ieta,ieps,ia]=c
                    val0[it,ieta,ieps,ia]=wvalue(a1)
                    
    return [lopt0,copt0,aopt0,val0]

#----------------------------------------------------------------------

#get steady value

#----------------------------------------------------------------------
def wagerate(x,y):
	w=(1-alpha)*x**alpha*y**(-alpha);
	return w


def interest(x,y):
	r=alpha*x**(alpha-1)*y**(1-alpha)-delta;
	return r

def getvaluess(kbar,lbar,averagehours,trbar0,taub,tauw0):
    
    global ia,ieta,ieps,it
    global asset0
    
    wbar=wagerate(kbar,lbar)
    rbar=interest(kbar,lbar)
    #ord = nage | neta | neps | na 
    #ord2 = t | neta | neps | na
    penbar=pen_repl*(1-taub-tauw0)*wbar*averagehours*np.ones((neps,1))
    
    wseq1=np.ones((nage,1))*wbar
    rseq1=np.ones((nage,1))*rbar
    penseq1=np.ones((nage,1))*penbar.T
    trseq1=np.ones((nage,1))*trbar0
    tauwseq1=np.ones((nage,1))*tauw0
    taubseq1=np.ones((nage,1))*taub
    
    ga=np.zeros((nage,neta,neps,na))
    
    #----------------------------
    
    #Age Path
    
    #----------------------------
    
    agen=np.zeros((nage,2))	
    cgen=np.zeros((nage,2))	
    lgen=np.zeros((t,1))
    lgeneps=np.zeros((t,2))	
    
    #----------------------------
    
    #Distribution
    
    #----------------------------
    
    fa=np.zeros((na,1))
    fnetincome=np.zeros((nin,1))
    fgrossincome=np.zeros((nin,1))
    fwage=np.zeros((nin,1))	
    fpension=np.zeros((nin,1))
    fconsumption=np.zeros((nin,1))
    
    #----------------------------
    
    #Economy Aggregate variable
    
    #----------------------------
    totalmeasure=0
    bigl=0
    bigcontrib=0
    bigc=0
    biga=0
    bigtax=0
    bigpensions=0
    
    [lopt,copt,aopt,v]=getoptpolicy(wseq1,rseq1,penseq1,trseq1,tauwseq1,taubseq1)
    
    #Initial State
    ga[0,0,0,na1]=1/4*mass[0]
    ga[0,1,0,na1]=1/4*mass[0]
    ga[0,0,1,na1]=1/4*mass[0]
    ga[0,1,1,na1]=1/4*mass[0]
    
    #--------------------------------------------------------------------
                    
    #From 1~T-1
                    
    #--------------------------------------------------------------------
                    
    for it in range(nage-1):
        for ia in range(na):
            asset0=agrid[ia]
            for ieps in range(neps):
                for ieta in range(eta):
                    measure=ga[it,ieta,ieps,ia]                    
                    
                    c=copt[it,ieta,ieps,ia] 
                    bigc=bigc+c*measure
                    
                    #-------------------------------------------------
                    
                    #change the statedy distribution of c
                    
                    #-------------------------------------------------
                    
                    if c<=0:
                        fconsumption[0]=fconsumption[0]+measure
                    elif c>=incomemax:
                        fconsumption[-1]=fconsumption[-1]+measure
                    else:
                        ini1=sum(c>incomegrid)
                        lam=(incomegrid[ini1]-c)/(incomegrid[ini1]-incomegrid[ini1-1])
                        fconsumption[ini1-1]=fconsumption[ini1-1]+lam*measure
                        fconsumption[ini1]=fconsumption[ini1]+(1-lam)*measure
                        
                    fa[ia]=fa[ia]+measure
                    
                    if it<t:
                        l0=lopt[it,ieta,ieps,ia]
                    
                    a1=aopt[it,ieta,ieps,ia]
                    
                    cgen[it,ieps]=cgen[it,ieps]+measure*c
                    agen[it,ieps]=agen[it,ieps]+measure*asset0
                    
                    bigtax=bigtax+taur*rbar*asset0*measure
                    totalmeasure=totalmeasure+measure
                    
                    
                    if it<t:
                        lgen[it]=lgen[it]+measure*l0
                        lgeneps[it,ieps]=lgeneps[it,ieps]+l0*measure
                        
                        bigl=bigl+l0*eta[ieta]*eps[ieps]*efage[it]*measure
                        bigtax=bigtax+tauw0*l0*efage[it]*eps[ieps]*eta[ieta]*wbar*measure
                        bigcontrib=bigcontrib+taub*l0*efage[it]*eps[ieps]*eta[ieta]*wbar*measure
                        
                        wageincome=l0*efage[it]*eps[ieps]*eta[ieta]*wbar
                        netincome=(1-tauw0-taub)*wageincome+trbar+(1-taur)*rbar*asset0
						
                        if wageincome<=0:
                            fwage[0]=fwage[0]+measure
                        elif wageincome>=incomemax:
                            fwage[-1]=fwage[-1]+measure
                        else:
                            ini1=sum(wageincome>incomegrid)
                            lam=(incomegrid[ini1]-wageincome)/(incomegrid[ini1]-incomegrid[ini1-1])
                            fwage[ini1-1]=fwage[ini1-1]+lam*measure
                            fwage[ini1]=fwage[ini1]+(1-lam)*measure
                            
                        grossincome=wageincome+trbar0+rbar*asset0
                        
                    else:
                        bigpensions=bigpensions+penbar[ieps]*measure
                        pensionincome=penbar[ieps]
                        grossincome=penbar[ieps]+trbar0+rbar*asset0
                        netincome=penbar[ieps]+trbar0+(1-taur)*rbar*asset0
                        
                    biga=biga+agrid[ia]*measure
                     
                    #-------------------------------------------------
                    
                    #change the statedy distribution of income
                    
                    #-------------------------------------------------
                            
                    if grossincome<=0:
                        fgrossincome[0]=fgrossincome[0]+measure
                    elif grossincome>=incomemax:
                        fgrossincome[nin-1]=fgrossincome[nin-1]+measure
                    else:
                        in1=sum(grossincome>incomegrid)
                        lam=(incomegrid[in1]-grossincome) / (incomegrid[in1]-incomegrid[in1-1] )
                        fgrossincome[in1-1]=fgrossincome[in1-1]+lam*measure
                        fgrossincome[in1]=fgrossincome[in1]+(1-lam)*measure
						
                    #-------------------------------------------------
                    
                    #change the statedy distribution of net income
                    
                    #-------------------------------------------------
                      
                    if netincome<=0:
                        fnetincome[0]=fnetincome[0]+measure
                    elif netincome>=incomemax:
                        fnetincome[nin-1]=fnetincome[nin-1]+measure
                    else:
                        in1=sum(netincome>incomegrid)
                        lam=(incomegrid[in1]-netincome) / (incomegrid[in1]-incomegrid[in1-1] )
                        fnetincome[in1-1]=fnetincome[in1-1]+lam*measure
                        fnetincome[in1]=fnetincome[in1]+(1-lam)*measure
                    
                    #-------------------------------------------------
                    
                    #change the statedy distribution of pension income
                    
                    #-------------------------------------------------
                    
                    
                    if it>=t:
                        if pensionincome<=0:
                            fpension[1]=fpension[1]+measure
                        elif pensionincome>=incomemax:
                            fpension[nin-1]=fpension[nin-1]+measure
                        else:
                            in1=sum(pensionincome>incomegrid)
                            lam=(incomegrid[in1]-pensionincome) / (incomegrid[in1]-incomegrid[in1-1] )
                            fpension[in1-1]=fpension[in1-1]+lam*measure
                            fpension[in1]=fpension[in1]+(1-lam)*measure
							
                    #-------------------------------------------------
                    
                    #change the statedy distribution of ga
                    
                    #-------------------------------------------------
                    
                    if a1<=assetmin:
                        ga[it+1,0,ieps,0]=ga[it+1,0,ieps,0]+pi_eta[ieta,0]*sp1[it]/(1+gn)*measure
                        ga[it+1,1,ieps,0]=ga[it+1,1,ieps,0]+pi_eta[ieta,1]*sp1[it]/(1+gn)*measure
					
                    elif a1>=assetmax:
                        ga[it+1,0,ieps,na-1]=ga[it+1,0,ieps,na-1]+pi_eta[ieta,0]*sp1[it]/(1+gn)*measure
                        ga[it+1,1,ieps,na-1]=ga[it+1,1,ieps,na-1]+pi_eta[ieta,1]*sp1[it]/(1+gn)*measure
					
                    else:
                        ia1=sum(agrid<a1)
                        lambda1=(agrid[ia1]-a1) / (agrid[ia1]-agrid[ia1-1] )
						
                        ga[it+1,0,ieps,ia1-1]=ga[it+1,0,ieps,ia1-1]	+lambda1 * pi_eta[ieta,0]*sp1[it]/(1+gn)*measure							
                        ga[it+1,1,ieps,ia1-1]=ga[it+1,1,ieps,ia1-1]	+lambda1 * pi_eta[ieta,1]*sp1[it]/(1+gn)*measure
                        ga[it+1,0,ieps,ia1]=ga[it+1,0,ieps,ia1]	+ (1-lambda1) * pi_eta[ieta,0]*sp1[it]/(1+gn)*measure
                        ga[it+1,1,ieps,ia1]=ga[it+1,1,ieps,ia1]+ (1-lambda1) * pi_eta[ieta,1]*sp1[it]/(1+gn)*measure

    #--------------------------------------------------------------------
                    
    #The Last Period
                    
    #--------------------------------------------------------------------
    
    it=nage-1
    
    for ia in range(na):
        asset0=agrid[ia]
        for ieps in range(neps):
            for ieta in range(neta):
                
                measure=ga[it,ieta,ieps,ia]
                totalmeasure=totalmeasure+measure
                 
                c=copt[it,ieta,ieps,ia]
                bigc=bigc+c*measure
                
                if c<=0:
                    fconsumption[0]=fconsumption[0]+measure
                elif c>=incomemax:
                    fconsumption[-1]=fconsumption[-1]+measure
                else:
                    ini1=sum(c>incomegrid)
                    lam=(incomegrid[ini1]-c)/(incomegrid[ini1]-incomegrid[ini1-1])
                    fconsumption[ini1-1]=fconsumption[ini1-1]+lam*measure
                    fconsumption[ini1]=fconsumption[ini1]+(1-lam)*measure
                
                fa[ia]=fa[ia]+measure
                
                grossincome=penbar[ieps]+trbar0+rbar*asset0
                netincome=penbar[ieps]+trbar0+(1-taur)*rbar*asset0

                pensionincome=penbar[ieps]                                                                                              				
                
                #-------------------------------------------------
                    
                #change the statedy distribution of income
                    
                #-------------------------------------------------
                            
                if grossincome<=0:
                    fgrossincome[0]=fgrossincome[0]+measure
                elif grossincome>=incomemax:
                    fgrossincome[nin-1]=fgrossincome[nin-1]+measure
                else:
                    in1=sum(grossincome>incomegrid)
                    lam=(incomegrid[in1]-grossincome) / (incomegrid[in1]-incomegrid[in1-1] )
                    fgrossincome[in1-1]=fgrossincome[in1-1]+lam*measure
                    fgrossincome[in1]=fgrossincome[in1]+(1-lam)*measure
					
                #-------------------------------------------------
                
                #change the statedy distribution of net income
                
                #-------------------------------------------------
                  
                if netincome<=0:
                    fnetincome[0]=fnetincome[0]+measure
                elif netincome>=incomemax:
                    fnetincome[nin-1]=fnetincome[nin-1]+measure
                else:
                    in1=sum(netincome>incomegrid)
                    lam=(incomegrid[in1]-netincome) / (incomegrid[in1]-incomegrid[in1-1] )
                    fnetincome[in1-1]=fnetincome[in1-1]+lam*measure
                    fnetincome[in1]=fnetincome[in1]+(1-lam)*measure
                
                #-------------------------------------------------
                
                #change the statedy distribution of pension income
                
                #-------------------------------------------------
                
                
                if it>=t:
                    if pensionincome<=0:
                        fpension[1]=fpension[1]+measure
                    elif pensionincome>=incomemax:
                        fpension[nin-1]=fpension[nin-1]+measure
                    else:
                        in1=sum(pensionincome>incomegrid)
                        lam=(incomegrid[in1]-pensionincome) / (incomegrid[in1]-incomegrid[in1-1] )
                        fpension[in1-1]=fpension[in1-1]+lam*measure
                        fpension[in1]=fpension[in1]+(1-lam)*measure
				
                cgen[it,ieps]=cgen[it,ieps]+measure*c
                agen[it,ieps]=agen[it,ieps]+measure*asset0
                bigtax=bigtax+taur*rbar*asset0*measure		# interest rate tax 75-year old
				
                bigpensions=bigpensions+penbar[ieps]*measure
                biga=biga+asset0*measure
                
    biga=sum(sum(agen))
    bigk=biga
    
    hoursnew=sum(lgen)/sum(mass)
    
    trnew=bigtax0+bigbequests-gbar
    tauwnew=tauw
    
    return [bigk,bigl,hoursnew,trnew,taubnew,tauwnew]
