import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as inter
plt.style.use('ggplot')
from scipy import linalg

#'''
def Slove_Quadratic_m2(A,B,C):
    def Slove_homogeneou_equation(X):
        def change(X):
            for j in range(np.shape(X)[0]):
                for i in range(j,np.shape(X)[0]):
                    #print(i,j)
                    if X[i][j]!=0:
                        X[[i,j],:]=X[[j,i],:]
                        #print(X)
            return X
        for x in range(np.shape(X)[0]-1):
            for j in range(np.shape(X)[1]-x-1):
                #print(np.shape(X)[1]-x-1)
                if X[x][x]==0:
                    continue
                X[j+x+1]=X[j+x+1]-X[x]*(X[j+x+1][x]/X[x][x])
        X=change(X)
        def Solve(X):
            H=np.zeros(np.shape(X)[0])
            if X[-1][-1]==0:
                H[-1]=1
            else:
                H[-1]=0
            for i in range(np.shape(X)[0]):
                #print('order:',np.shape(X)[0]-i-1)
                if X[np.shape(X)[0]-i-1][np.shape(X)[0]-i-1]==0:
                    H[np.shape(X)[0]-i-1]=1
                    #print('H:',H)
                else:
                    #print('X',X[np.shape(X)[0]-i-1][np.shape(X)[0]-i:])
                    #print('H',H[np.shape(X)[0]-i:])
                    S=-np.array(X[np.shape(X)[0]-i-1][np.shape(X)[0]-i:])@np.array(H[np.shape(X)[0]-i:])
                    H[np.shape(X)[0]-i-1]=S/X[np.shape(X)[0]-i-1][np.shape(X)[0]-i-1]
            return H
        L=Solve(X)
        #print('Zero',np.dot(X,L))
        #print(X)
        return L
    def find_root(X,Y):
        plt.plot(X,Y)
        slove=[]
        for y in Y:
            if abs(y)<1e-3:
                
                if len(slove)==0:
                    slove.append(X[Y.index(y)])
                else:
                    if abs(X[Y.index(y)]-slove[-1])>1:
                        slove.append(X[Y.index(y)])
        return slove

    
    lambd=np.linspace(-1,1,25000)
    def f_lambd(lambd,A,B,C):
        Y=[]
        for i in lambd:
            Y.append(np.linalg.det(A*i**2+B*i+C))
        return Y
    
    Y=f_lambd(lambd,A,B,C)
    Eigenvalue=find_root(lambd,Y)
    print(Eigenvalue)
    
    def get_vector(Eigenvalue,A,B,C):
        S=[]
        for i in Eigenvalue:
            X=A*i**2+B*i+C
            X=Slove_homogeneou_equation(X)
            S.append(X)
            #print(X)
            
            #print(np.dot(X,orth(X).T[1]))
        S=np.array(S)
        return S
    
    Orth=get_vector(Eigenvalue,A,B,C).T
    #print(Orth)
    
    #P=np.dot(Orth,np.dot(np.diag(Eigenvalue),np.linalg.inv(Orth)))
    P=Orth@Eigenvalue@(1/Orth)
    print(P)
    #print(A@P@P+B@P+C)
    return P

def Slove_Quadratic(A,B,C):
    Monit_var=0
    def Dispose_DE(A,B,C):
        UP=np.column_stack((B,C))
        n=np.shape(B)[0]
        I=np.identity(n)
        #print(I)
        O=np.zeros((n,n))
        #print(O)
        DOWN=np.column_stack((I,O))
        #print(S)
        D=np.row_stack((UP,DOWN))
        UP=np.column_stack((A,O))
        DOWN=np.column_stack((O,I))
        E=np.row_stack((UP,DOWN))
        return D,E
    D,E=Dispose_DE(A,B,C)
    
    #Eigenvalue and EigenVector:Generalized Right Eigenvalue Problem
    X=linalg.eig(D,E)
    
    Monit_var=X
    Comprise1=D@X[1]
    Comprise2=E@X[1]@np.diag(X[0])
    #print(X[0],'\n',X[1])
    def arrange_the_vector(M):
        
        sequence=np.argsort(abs(M[0]))
        n=len(M[0])
        se=np.ones(n)
        M_a=np.zeros(shape=(n,len(M[1][0])))
        M1=np.array(M[1])
        M1=M1.T
        for i in range(n):
            M_a[i]=M1[sequence[i]]
            se[i]=M[0][sequence[i]]
        M_a=M_a.T
        return M_a,se
    def arrange_the_vec(M):
        n=len(M[0])
        se=np.ones(n)
        number=0
        M_a=np.zeros(shape=(n,len(M[1][0])))
        M1=np.array(M[1])
        M1=M1.T
        for i in range(n):
            if abs(M[0][i])<1:
                se[number]=M[0][i]
                M_a[number]=M1[i]
                number=number+1
        for i in range(n):
            if abs(M[0][i])>1:
                se[number]=M[0][i]
                M_a[number]=M1[i]
                number=number+1
        M_a=M_a.T
        return M_a,se
    X,sequence=arrange_the_vector(X)
    X_comp=X
    sequence_comp=sequence
    #print('\n',sequence,'\n',X)
    
    def M_split(X):
        X=np.array(X)
        X_UP,X_DOWN=np.hsplit(X,2)
        X11,X12=np.vsplit(X_UP,2)
        X21,X22=np.vsplit(X_DOWN,2)
        return X11,X21,X12,X22
    #Delta=X^-1*E^-1*D*X
    X11,X21,X12,X22=M_split(X)
    #print(X11,'\n',sequence[0:3]*X12)
    #Delta_one=np.dot(np.linalg.inv(X12),X11)
    Delta_one=np.diag(sequence[0:np.shape(A)[0]])
    
    #P=np.dot(X12,np.dot(Delta_one,np.linalg.inv(X12)))
    #print(X11,X12)
    P=X11@np.linalg.inv(X12)
    return P,Monit_var,D,E,Comprise1,Comprise2,X_comp,sequence_comp,X11,X12
#'''  
'''
def Slove_Quadratic_m3(A,B,C):
    P=np.ones((np.shape(A)))
    n=200
    
    for i in range(n):
        
    return P
'''
#Kronecker Product
def Kronecker(A,B,A_T=0,B_T=0):
    if len(A.shape)==1:
        if A_T==1:
            A=np.vstack(A)
        else:
            A=np.vstack(A).T
    #print('As',A.shape)
    if len(B.shape)==1:
        if B_T==1:
            B=np.vstack(B)
        else:
            B=np.vstack(B).T
    #print(A,B)  
    lengthA=A.shape[0]
    lengthB=B.shape[0]
    widthA=A.shape[1]
    widthB=B.shape[1]
    
    C=np.zeros((lengthA*lengthB,widthA*widthB))
    #print(C.shape)
    #print('lenA',lengthA,'lenB',lengthB)
    for i in range(lengthA):
        for j in range(widthA):
            for ib in range(lengthB):
                for jb in range(widthB):
                    #print(i*length+ib,j*widthA+jb)
                    C[i*lengthB+ib,j*widthA+jb]=A[i,j]*B[ib,jb]
    return C

#Dispose the matrix
def vec(A):
    if len(np.shape(A))==1:
        A=np.vstack(A).T
    V=[]
    for i in range(np.shape(A)[0]):
        for j in range(np.shape(A)[1]):
            V.append(A[i][j])
    return np.array(V).T

def Slove_DSGE(A,B,C,D,F,J,K,L,G,H,M,N):
    #x_t=Px_{t-1}+Qz_t
    #y_t=Rx_{t-1}+Sz_t
    #P is an variable
    #print(H-np.dot(K,np.dot(np.linalg.inv(C),B)))
    if len(H-np.dot(K,np.dot(np.linalg.inv(C),B)))==1:
        B_sq=-(np.dot(J,np.dot(np.linalg.inv(C),B))-G+np.dot(K,np.dot(np.linalg.inv(C),A)))
        #print(B_sq**2)
        #print(np.dot(B_sq,B_sq))
        A_sq=(F-np.dot(J,np.dot(np.linalg.inv(C),A)))
        C_sq=H-np.dot(K,np.dot(np.linalg.inv(C),B))
        P_pos=(-B_sq+np.sqrt(np.dot(B_sq,B_sq)-4*np.dot(A_sq,C_sq)))/2/A_sq
        P_neg=(-B_sq-np.sqrt(np.dot(B_sq,B_sq)-4*np.dot(A_sq,C_sq)))/2/A_sq
        P=[P_pos,P_neg]
        print('P:',P)
    else:
        P,Monit_var=Slove_Quadratic(A,B,C)
        
        #print('P:',P)
    #P=min(P)
    #==================================================#
    #R
    #print(P,np.dot(A,2))
    R=-np.dot(np.linalg.inv(C),(np.dot(A,P[0])+B))
    #print(R)
    #==================================================#
    #Q
    I_k=np.ones(N.shape[0])
    
    K_A=(Kronecker(N.T,F-np.dot(np.dot(J,np.linalg.inv(C)),A),1,0)+
         Kronecker(I_k,np.dot(J,R)+np.dot(F,P)+
                   G-np.dot(K,np.dot(np.linalg.inv(C),A))))
    Q=np.linalg.inv(K_A)*vec(np.dot((np.dot(J,np.dot(np.linalg.inv(C),D))-L),N)+ \
                    np.dot(K,np.dot(np.linalg.inv(C),D))-M)
    
    #==================================================#
    S=-np.dot(np.linalg.inv(C),(np.dot(A,Q[0][0])+D))
    return [P,Q[0],R,S]
#'''
r_bar=0.0351
Y_bar=1.1875
H_bar=0.3206
K_bar=12.1795
C_bar=0.8830

beta=0.99
delta=0.025
theta=0.36
A=1.72
rho_w=0.7
phi_w=21
gamma=0.95
pi=0.48

A=np.array([
            [K_bar,C_bar,-C_bar,0],
            [0,0,1,0],
            [0,0,0,0],
            [0,0,0,0]
            ])
B=np.array([
            [-(1-delta)*K_bar,0,0,0],
            [0,-1,0,0],
            [1,0,0,0],
            [-theta,0,0,0],
            ])
C=np.array([
            [0,0,-Y_bar,0],
            [0,1,0,0],
            [1,0,-1,0],
            [0,0,1,-(1-theta)]
            ])
D=np.array([
            [0,0],
            [0,-1],
            [0,0],
            [-1,0]
            ])
F=np.array([
            [0,1,-2,0],
            [0,0,0,0],
            [0,0,-(1-rho_w)*(1-beta*rho_w),-beta*rho_w],
            [0,0,0,0]
            ])
G=np.array([[0,0,1,0],
            [0,1,0,0],
            [0,0,0,1+beta*rho_w**2+(1-beta*rho_w)*phi_w*rho_w*H_bar/(1-H_bar)],
            [0,0,-1,1]
            ])
H=np.array([[0,0,0,0],
            [0,-1,0,0],
            [0,0,0,-rho_w-(1-beta*rho_w)*phi_w*rho_w*H_bar/(1-H_bar)],
            [0,0,0,0]])
J=np.array([[-beta*r_bar,-1,0,0],
            [0,0,0,0],
            [0,-(1-rho_w)*(1-beta*rho_w),0,0],
            [0,0,0,0]])

K=np.array([
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,-(1-rho_w)*(1-beta*rho_w)*H_bar/(1-H_bar)],
            [0,0,-1,1]
            ])
L=np.array([[0,pi],
            [0,0],
            [0,0],
            [0,0]])
M=np.array([[0,0],
            [0,-1],
            [0,0],
            [0,0]
            ])
N=np.array([[gamma,0],
            [0,pi]])
#print(len(D.shape))
P,Monit_var,D,E,Comprise1,Comprise2,X_comp,sequence_comp,X11,X12=Slove_Quadratic(A,B,C)
#P=Slove_Quadratic_m3(A,B,C)
#Slove=Slove_DSGE(A=A,B=B,C=C,D=D,F=F,G=G,H=H,J=J,K=K,L=L,M=M,N=N)
'''
#'''
Times=28
K_t=[0]
lambda_t=np.zeros(Times)
lambda_t[2]=lambda_t[2]+0.01
Y_t=[]
C_t=[]
H_t=[]
r_t=[]
for i in range(Times):
    K_t.append(Slove[0][0]*K_t[-1]+Slove[1][0]*lambda_t[i])
for i in range(Times):
    Y_t.append(Slove[2][0]*K_t[i]+Slove[3][0]*lambda_t[i])
    C_t.append(Slove[2][1]*K_t[i]+Slove[3][1]*lambda_t[i])
    H_t.append(Slove[2][2]*K_t[i]+Slove[3][2]*lambda_t[i])
    r_t.append(Slove[2][3]*K_t[i]+Slove[3][3]*lambda_t[i])
X_grid=np.linspace(0,Times-1,300)
plt.subplot(321)
K_t=inter.interp1d(range(0,Times+1),K_t,kind='quadratic')
plt.plot(X_grid,K_t(X_grid))
plt.title('$K$')

plt.subplot(322)
Y_t=inter.interp1d(range(0,Times),Y_t,kind='quadratic')
plt.plot(X_grid,Y_t(X_grid))
plt.title('$Yield$')

plt.subplot(323)
C_t=inter.interp1d(range(0,Times),C_t,kind='quadratic')
plt.plot(X_grid,C_t(X_grid))
plt.title('$Consumption$')

plt.subplot(324)
H_t=inter.interp1d(range(0,Times),H_t,kind='quadratic')
plt.plot(X_grid,H_t(X_grid))
plt.title('$Human$')

plt.subplot(325)
r_t=inter.interp1d(range(0,Times),r_t,kind='quadratic')
plt.plot(X_grid,r_t(X_grid))
plt.title('$r$')

#'''

