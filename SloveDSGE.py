import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as inter

plt.style.use('ggplot')

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
#'''
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
    B_sq=-(np.dot(J,np.dot(np.linalg.inv(C),B))-G+np.dot(K,np.dot(np.linalg.inv(C),A)))
    A_sq=(F-np.dot(J,np.dot(np.linalg.inv(C),A)))
    C_sq=H-np.dot(K,np.dot(np.linalg.inv(C),B))
    P_pos=(-B_sq+np.sqrt(B_sq**2-4*np.dot(A_sq,C_sq)))/2/A_sq
    P_neg=(-B_sq-np.sqrt(B_sq**2-4*np.dot(A_sq,C_sq)))/2/A_sq
    P=[P_pos,P_neg]
    #print(P)
    P=min(P)
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
A=np.array([0,-12.6698,0,0]).T
A_T=1
B=np.array([0,12.3530,0.36,-1]).T
B_T=1
C=np.array([[1,-1,-1.5004,0],
            [1.2353,-0.9186,0,0],
            [-1,0,0.64,0],
            [1,0,0,-1]])
D=np.array([0,0,1,0]).T
D_T=1
F=np.array([0])
F_T=0
G=np.array([0])
G_T=0
H=np.array([0])
H_T=0
J=np.array([0,-1,0,0.0348])
J_T=0
K=np.array([0,1,0,0])
K_T=0
L=np.array([0])
L_T=0
M=np.array([0])
M_T=0
N=np.array([0.95])
N_T=0
#print(len(D.shape))
Slove=Slove_DSGE(A=A,B=B,C=C,D=D,F=F,G=G,H=H,J=J,K=K,L=L,M=M,N=N)
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







