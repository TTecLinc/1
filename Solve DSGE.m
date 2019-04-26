%The E
% theta=0.36;
% beta=0.99;
% gamma=0.95;
% delta=0.025;
% pie=0.48;
% varphi=1.3180;
% eta=0.8;
% BB=-2.5805;
% rbar=1/beta-(1-delta);
% wbar=(1-theta)*(theta/rbar)^(theta/(1-theta));
% gbar=(-wbar*beta/(BB*varphi))^(1/eta)*(varphi-1);
% KoverH=(theta/rbar)^(1/(1-theta));
% pbar=(-BB/(wbar*beta))^(1/eta)*(varphi)^(1/eta-1);
% hibar=1/(pbar*(wbar+(rbar-delta)*KoverH));
% Hbar=hibar;
% kibar=KoverH*Hbar;
% Kbar=kibar;
% cbar=1/(pbar*varphi);
% mbar=1;
% %build the matrices 
% A=[0 -1/varphi
%     Kbar 0
%     0 0
%     0 0];
% B=[0 0
%     -(rbar+1-delta)*Kbar 0
%     1-theta 0
%     -theta 0];
% C=[0 0 pbar*gbar 0
%     -rbar*Kbar -wbar*Hbar -1/pbar -wbar*Hbar
%     1 0 0 -(1-theta)
%     0 1 0 theta];
% D=[0 pbar*gbar
%     0 0
%     -1 0
%     -1 0];
% F=[0 0
%     0 (eta-1)];
% G=[0 0
%     0 0];
% H=[0 0
%     0 0];
% J=[beta*rbar -1 0 0
%     0 0 (eta-1) 0];
% K=[0 1 0 0
%     0 1 1 0];
% L=[0 0
%     0 0];
% M=[0 0
%     0 0];
% I1=[1 0
%     0 1];
% Z1=zeros(2);
% N=[gamma 0
%     0 pie];
theta=0.36;
beta=0.99;
gamma=0.95;
delta=0.025;
pie=0.48;
rho_w=0.7;
psi_w=21;
r_bar=0.0351
H_bar=0.3206;
K_bar=12.1795;
Y_bar=1.1875;
C_bar=0.8830;
%build the matrices 
A=[K_bar C_bar -C_bar 0
    0 0 1 0
    0 0 0 0
    0 0 0 0];
B=[-(1-delta)*K_bar 0 0 0
    0 -1 0 0
    1 0 0 0
    -theta 0 0 0];
C=[0 0 -Y_bar 0
    0 1 0 0
    1 0 -1 0
    0 0 1 -(1-theta)];
D=[0 0
    0 -1
    0 0
    -1 0];
F=[0 1 -2 0
    0 0 0 0
    0 0 -(1-rho_w)*(1-beta*rho_w) -beta*rho_w
    0 0 0 0];
G=[0 0 1 0
    0 1 0 0
    0 0 0 1+beta*rho_w^2+(1-beta*rho_w)*psi_w*rho_w*H_bar/(1-H_bar)
    0 0 -1 1];
H=[0 0 0 0
    0 -1 0 0
    0 0 0 -rho_w-(1-beta*rho_w)*psi_w*rho_w*H_bar/(1-H_bar)
    0 0 0 0];
J=[-r_bar*beta -1 0 0
    0 0 0 0
    0 -(1-rho_w)*(1-beta*rho_w) 0 0
    0 0 0 0];
K=[0 0 0 0
    0 0 0 0
    0 0 0 -(1-rho_w)*(1-beta*rho_w)*H_bar/(1-H_bar)
    0 0 -1 1];
L=[0 pie
    0 0
    0 0
    0 0];
M=[0 0
    0 -1
    0 0
    0 0];
I1=eye(size(J));
Z1=zeros(size(J));
N=[gamma 0
    0 pie];
%=========================================================================
%Find the solution to the matrix quadratic 
invC=inv(C);
psy=F-J*invC*A;
lambda=J*invC*B-G+K*invC*A;
T=K*invC*B-H;
AA1=[lambda T
    I1 Z1];
AA2=[psy Z1
    Z1 I1];
[eigvec,eigval]=eig(AA1,AA2);
diageigval=diag(eigval);
%Select the stable egivenvalue
iz=find(abs(diageigval)<1.05);
max_num_vec=size(J,2)
%iz=iz(1:max_num_vec);
DD=diag(diageigval(iz));
ei=size(eigvec);
EE=eigvec(ei/2+1:ei,iz);
P=EE*DD*inv(EE)
%===============================================================%
R=-invC*(A*P+B)
I2=[1 0
    0 1];
QQ=kron(N,(F-J*invC*A))+kron(I2,(J*R+F*P+G-K*invC*A));
invQQ=inv(QQ);
QQQ=((J*invC*D-L)*N+K*invC*D-M);
[aa,bb]=size(QQQ);
Qfindvert=[];
for ij=1:bb
    Qfindvert=vertcat(Qfindvert,QQQ(:,ij));
end
Qvert=invQQ*Qfindvert;
Q=[];
for ij=1:bb
    begini=(ij-1)*aa+1;
    endi=ij*aa;
    Q=[Q Qvert(begini:endi,1)];
end
S=-invC*(A*Q+D);







