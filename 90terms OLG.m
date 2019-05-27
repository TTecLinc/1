%   Ch6_optimal_pension.m
%   computes the transition dynamics in Ch. 6.4 of Heer (2018)
%   VALUE FUNCTION ITERATION 
%
%   author: Burkhard Heer
%   date: July 7, 2018
%
%	Algorithm: See Algorithm 9.2.1 in Heer/Maußner DGE Modeling
%
%	ASSUMPTION: The policy change is implemented gradually 
%	(linearly between the years year_change and year_policy_final)
%
%
%	First run the program for case_ss=1 and _policy_change=0 
%	to save steady state values and benchmark transition (which saves lifetime utility of transition generations)
%
%	Next run the program for case_ss=0 and _policy_change=1 

clear all     % clear variables and functions from memory
close all     % closes all figure windows
clc           % clears the command window and homes the cursor.


disp('This program computes the solution to the Optimal Pension Problem');
disp(' ' );
disp('in the Chapter 6.4 of Heer, Public Economics - The Macroeconomic Perspective.');
disp(' ' );
disp('Hit any key when ready.....');
% pause;
tic;        % computational time

% define global variables
def_global

%
%   STEP 1: Parameterization and calibration
%
% Parameterization
fhandle_function1 = @value1;
fhandle_function2 = @wvalue;
taub=0.05;
pen_repl=0.5;       %  pension net replacement rate
pen_repl_old=0.5;
pen_repl_new=0.5;
% pen_repl_new=0.20;
% pen_repl_new=0.33;
beta1=1.011;
r=0.032;            % initial value of the interest rate */
sigma=2;            % coefficient of relative risk aversion */
alpha1=0.36;         % production elasticity of capital */
delta=0.08;         % rate of depreciation */
gam=0.31;           % utility rshare of consumption 
gy=0.195;           % government consumption share
tauw=0.248;         % labor income tax rate
taur=0.429;         % interest income tax rate
garate=0.02;		% productivity growth rate
% garate=0;         % no growth
laborexogenous=0.3;
tr=25;              % retirement periods
t=45;               % working time 
nage=t+tr;
neta = 2;
neps = 2;
% calibration of the efficiency types
eps1=[0.57; 1.43];
eta1=[0.727; 1.273];
pi_eta=[0.98 0.02;0.02 0.98];

year_change=2015;       % year when the policy change is announced
year_change1=2015;      % year when policy change starts
year_final=2100;        % year when the population is assumed to be stationary
year_finalss=2050;      % computation of steady state for comparative steady-state analysis
calib=1;               % 1 - during first run, to calibrate gbar 
case_ss=0;              % 1 -- computation of the steady states (only during first run), 0 -- no computation
% ssfinal=1;             % do not change, just for checking
policy_change=1;       % 0 -- benchmark, 1 -- policy change 
case_endogenouslabor=1; % 1 - endogenous labor, 2 - exogenous labor
case_pension=2;		% 1 -- proportional to permanent type
						% 2 -- lump-sum
case_UNscen=1;          % 1 -- medium variant, 2 -- low variant, 3 -- high variant for population projection UN (2015)
endog_tr=1;             % 1 -- tranfers adjust to balance fiscal budget, 0 -- wage income tax adjusts
nage1=15;
year1=year_change;		% 2010 or 2050: population parameters, multiple of 5
clc;


% computational parameters
tol_golden = 1e-5;      % tolerance for golden section search
ntrans=250;             % number of total transition periods 
ntrans_finalss=(year_finalss-2015)+1;		% number of years between initial and final steady state
											% for comparative steady-state analysis
year_policy_final=year_change1+45;          % year when the whole policy change is implemented
                                            % workers have a lifetime to adjust
lmin=0.01;              % small constant for labor supply, to rule out corner solution in the worker's 
                        % optimization problem
lmax=0.90;              % maximum hours
neg=-1e10;              % initial value for value function */
psi1=0.001;              % small constant in utility function (to allow for zero consumption as the lower bound)
VI_method='linear';           % 'linear' --- linear interpolation, 'cubic' --- cubic spline 
% Vq = interp2(,method) X,Y,V,Xq,Yq) returns interpolated values of a function of two 
% variables at specific query points using linear/cubic interpolation. 
% The results always pass through the original sampling of the function. 
% X and Y contain the coordinates of the sample points. 
% V contains the corresponding function values at each sample point. 
% Xq and Yq contain the coordinates of the query points.
maxitrans=20;               % maximum number of iteration over transition pathes
maxitrans=1;
phi_updatetrans=0.8;        % updating parameter for vector of transition pathes
phig_updatetrans=0.95;      % updating parameter for vector of transition pathes
toltrans=0.001;             % maximum divergence of final path for K 
itrans=0;
krittrans=1+toltrans;

% 
% asset grid
assetmax=10;    % upper limit of asset grid on a 
na=50;          % number of grid points on assets
assetmin=0;
agrid=linspace(0,assetmax,na);   % asset grid	
agrid = agrid';     % column vector
% income grid
nin=1000;
incomemax=4;
incomegrid=linspace(0,incomemax,nin);   % income grid
incomegrid=incomegrid';     % column vector
cgen0=zeros(nage,neps);     % generational consumption

%
%   STEP 2: Preparation of the population data
%
% loading of the demographic data from survival_probs_US.xls";
survivalprobs=xlsread('survival_probs_US.xls',2,'F22:AI37');
survivalprobs=flipud(survivalprobs);
popgrowth=xlsread('survival_probs_US.xls',2,'F41:AI43');
timespan=linspace(1950,2095,30);
nrate=[timespan' popgrowth'];

% computation of movering averages of population data
movavperiods=4;         %  number of moving average periods for computation of population growth and suvival probs
year0=(year1-1950)/5+1;
if year0==round(year0)
	popgrowth=nrate(year0-movavperiods+1:year0,case_UNscen+1);
	popgrowth=mean(popgrowth);
	popgrowth=popgrowth/100;        % data is expressed in percentage numbers
	sp=survivalprobs(1:nage1,year0-movavperiods+1:year0);
	sp=mean(sp,2);
else
	disp('year1 must be a multiple of 5');
	pause;
end

year1
popgrowth
gn=popgrowth;

% conversion of 5-year rates to annual rates
% cubic spline interpolation
%
% vq = interp1(x,v,xq,'spline') returns interpolated values of a 1-D function at 
% specific query points using linear interpolation. 
% Vector x contains the sample points, and v contains the corresponding values, v(x). 
% Vector xq contains the coordinates of the query points.
%
age=linspace(20,20+5*(nage1-1),nage1);
nage0=5*(nage1-1)+1;
age1=linspace(20,20+5*(nage1-1),5*(nage1-1)+1);
age=age';       % column vector
age1=age1';     
sp1=interp1(age,sp,age1,'spline');
sp1=sp1.^(1/5);          % conversion from 5-year survival probs to 1-year probs

figure
xlabel('age');
title('5-annual survival probability');
plot(age,sp);
%pause;
figure
plot(age1,sp1);
title('Annual survival probability');
%pause;	

sp1=sp1(1:nage);
sp1initial=sp1;
gninitial=gn;

% loading of the age-efficiency profile, as estimated by Hansen
efage=xlsread('efficiency_profile.xls',1,'A1:A45');
efage=efage/mean(efage);
age2=linspace(20,64,45);
figure
plot(age2,efage);
xlabel('age');
title('Productivity');



% computation of cohort mass in year_change 
mass=ones(nage,1);

for i=2:1:nage
    mass(i)=mass(i-1)*sp1(i-1)/(1+gn);
end
mass=mass/sum(mass);
massinitial=mass;

figure
plot(age1(1:nage),mass);
xlabel('age');
title('Cohort mass');	

%
% Step 2: preparation of the annual survival rates and population growth rates
%


year1=year_change;				% 2015 or 2050: population parameters, multiple of 5
year0=(year1-1950)/5+1;
sptotal=sp1initial(1:nage);
gntotal=gninitial;

while year0<size(nrate,1)
	year0=year0+1;
	if year0==round(year0)
		popgrowth=nrate(year0-movavperiods+1:year0,case_UNscen+1);
		popgrowth=mean(popgrowth,1);
		popgrowth=popgrowth/100;	% data is expressed in percentage numbers
		sp=survivalprobs(1:nage1,year0-movavperiods+1:year0);
		sp=mean(sp,2);
    else
		disp('year1 must be a multiple of 5');
		pause;
    end
	year0
	popgrowth

	% conversion of 5-year rates to annual rates
	% cubic spline interpolation
    sp1=interp1(age,sp,age1,'spline');
    sp1=sp1.^(1/5);          % conversion from 5-year survival probs to 1-year probs
   
	sptotal=[sptotal sp1(1:nage)];
	gntotal=[gntotal; popgrowth];
end

%
% linear interpolation of population growth rates and survival probs
%

nsp=size(sptotal,2);
yeart=linspace(year1,year1+(nsp-1)*5,nsp);
rowsspall=(nsp-1)*5+1;
yeart1=linspace(year1,year1+rowsspall-1,rowsspall);

gnall=zeros(rowsspall,1);
spall=zeros(nage,rowsspall);

spall(1:nage,rowsspall)=sptotal(1:nage,nsp);
gnall(rowsspall)=gntotal(nsp);

i=0;
for i=1:1:nsp-1    
	gnall((i-1)*5+1)=gntotal(i);
	gnall((i-1)*5+2)=4/5*gntotal(i)+1/5*gntotal(i+1);
	gnall((i-1)*5+3)=3/5*gntotal(i)+2/5*gntotal(i+1);
	gnall((i-1)*5+4)=2/5*gntotal(i)+3/5*gntotal(i+1);
	gnall((i-1)*5+5)=1/5*gntotal(i)+4/5*gntotal(i+1);
	
	spall(1:nage,(i-1)*5+1)=sptotal(1:nage,i);
	spall(1:nage,(i-1)*5+2)=4/5*sptotal(1:nage,i)+1/5*sptotal(1:nage,i+1);
	spall(1:nage,(i-1)*5+3)=3/5*sptotal(1:nage,i)+2/5*sptotal(1:nage,i+1);
	spall(1:nage,(i-1)*5+4)=2/5*sptotal(1:nage,i)+3/5*sptotal(1:nage,i+1);
	spall(1:nage,(i-1)*5+5)=1/5*sptotal(1:nage,i)+4/5*sptotal(1:nage,i+1);
	
end

%
%		computation of the mass during the transition
%

massvec=zeros(nage,ntrans);			% mass of the generation j in period tp
massvec0=massvec;					% just to check if the mass is equal to one during each transition period
massvec1=zeros(ntrans,1);
mass0=massinitial;
massvec(1:nage,1)=mass0;
mass1=zeros(nage,1);
for i=1:1:ntrans-1	
	j=1;
	mass1(1)=mass0(1)*(1+gnproc(i+1));
    for j=2:1:nage
		mass1(j)=mass0(j-1)*surviveprob(j-1,i);
    end
	mass0=mass1;
	massvec(1:nage,i+1)=mass1;
end

massyear=sum(massvec,2);

disp('computation of massvec and spvec complete');
disp('dependency ratios');
dp=zeros(ntrans,1);
for i=1:1:ntrans
	dp(i)=sum(massvec(t+1:nage,i))/sum(massvec(1:t,i));
end

figure
plot(linspace(2015,2015+ntrans-1,ntrans),dp);
title('Dependency ratio 65+/(20-65)');
xlabel('Year');
% pause;




%
%   STEP 3: Computation of the initial steady state in year "year_change"
%
%


sp1=sp1initial;
gn=gninitial;

% initialization of the factor price sequences 
%	for the cohort optimization problem
%
wseq=zeros(nage,1);
rseq=zeros(nage,1);
trseq=zeros(nage,1);
penseq=zeros(nage,neps);
taurseq=zeros(nage,1);
tauwseq=zeros(nage,1);
taubseq=zeros(nage,1);
wuseq=zeros(nage,neps);
hickstrseq=zeros(nage,1);

aoptseq = zeros(nage,neps,neta,na);
vseq = zeros(nage,neps,neta,na);
coptseq = zeros(nage,neps,neta,na);
loptseq = zeros(t,na,neps,neta,na);
% ord2 = t | neta | neps | na ;
aopt=aoptseq;
lopt=loptseq;
copt=coptseq;
v=vseq;
ve=v; vu=v; 
% case_ui=0; ui_rate=0; ord1=0; ord3=0; zeta1=0;


%
% initialization as global variable
%
it=0;	
vr1=0;
vr11=0;
v11e=0;
v12e=0;
ieta=0; 
ieta1=0; 
ia=0; 
ia1=0;
asset0=0; 
ieps=0;
ga=zeros(nage,neps,neta,na);
% age-profiles
agen=zeros(nage,2);     % assets workers
cgen=zeros(nage,2);		% consumption epsilon=1,2
lgen=zeros(t,1);		% working hours workers
lgeneps=zeros(t,2);		% working hours, for each epsilon
% distribution functions for Gini coefficients
fa=zeros(na,1);         % distribution of wealth a
fnetincome=zeros(nin,1);    % net income
fgrossincome=zeros(nin,1);	% gross income
fwage=zeros(nin,1);			% wage income
fpension=zeros(nin,1);		%  pension income
fconsumption=zeros(nin,1);	% consumption
giniwage=0;
ginigrossincome=0;
gininetincome=0;
giniwealth=0;
ginipension=0;
giniconsumption=0;
progressivityindex=0;
bigc=0;
bigy=0;
welfare=0;
totalpension=0;


if case_ss==1
    % initial guess for l and k
    averagehours = 0.3;
    lbar = averagehours * sum(mass(1:t));
    kbar = (alpha1 /(1/beta1-1+delta) )^(1/(1-alpha1))*lbar;
    ybar = kbar^(alpha1)*lbar^(1-alpha1);
    gbar = gy * ybar;
    trbar = 0;
    taub = 0.1;
    close all;
    averagehours = 0.2980;
    lbar = 0.27463;
    kbar = 1.10991;
    ybar = kbar^(alpha1)*lbar^(1-alpha1);
    gbar = gy * ybar;
    trbar = 0.02293;
    taub = 0.08532;
    x0=[kbar lbar averagehours trbar taub tauw];		
    y=getss(x0); 
    y
ttime=toc;
disp('Computation of transition complete');
disp(['Elapsed time: ',num2str(ttime/60),' minutes']);
pause;
  
   

    xagg = fsolve(@(x)getss(x),x0);
    xagg = x0;
    clc;
    xagg
   
ttime=toc;
disp('Computation of transition complete');
disp(['Elapsed time: ',num2str(ttime/60),' minutes']);
pause;
  
    xinitial=xagg;
    kbarinitial=xagg(1);
    lbarinitial=xagg(2);
    averagehoursinitial=xagg(3);
    trbarinitial=xagg(4);
    taubinitial=xagg(5);
    tauwinitial=xagg(6);
    gainitial=ga;	% initial distribution of assets in period 0
    gainitial0=ga;
    valueinitial=v;
    
    
    massinitial=mass;
    gbar2015=gbar; 
    gbarinitial=gbar;
    loptinitial=lopt; 
    coptinitial=copt; 
    aoptinitial=aopt;
    % save steady state values to file, for future runs 
    save('Ch6_output.mat','kbarinitial', 'lbarinitial', 'averagehoursinitial', 'trbarinitial','taubinitial','tauwinitial');
    save('Ch6_output.mat','sp1initial', 'gninitial', 'massinitial', 'gbar2015', 'gainitial', '-append');
    save('Ch6_output.mat','loptinitial','coptinitial','aoptinitial','valueinitial', '-append');

    cgen01=cgen(:,1)./mass*2;
    cgen02=cgen(:,2)./mass*2;
    cgen0=[cgen01 cgen02];
    lgen0=lgen./mass(1:45);
    agen01=agen(:,1)./mass*2;
    agen02=agen(:,2)./mass*2;
    agen0=[agen01 agen02];
    lgeneps01=lgeneps(:,1)./mass(1:45)*2;
    lgeneps02=lgeneps(:,2)./mass(1:45)*2;
    lgeneps0 = [ lgeneps01 lgeneps02];

	periods=linspace(20,64,45);
	figure
	plot(periods,lgen0);
	xlabel('Age');
	ylabel('Average working hours');
%	pause;

    figure
    plot(periods,lgeneps0);
    legend('low productivity','high productivity');
	xlabel('Age');
	ylabel('Average working hours');
%	pause;
	
	disp('value newborn');
	v(1,1,1:2,1:2)
 %   pause;



	periods=linspace(20,20+nage-1,nage);
	
    figure
    plot(periods,agen0);
	ylabel('Average assets');
	xlabel('Age');
%	pause;
    
    figure
    plot(periods,cgen0);
	ylabel('Average consumption');
	xlabel('Age');
%	pause;
	
	

	disp('credit-constrained: '); 
    fa(1) 
 %   pause;


	if fa(1)<1
        figure
        plot(agrid,fa)
		xlabel('asset');
		ylabel('distribution');
%		pause;
    end

    figure 
    plot(incomegrid,fgrossincome);
	xlabel('gross income');
	ylabel('distribution');
	
%	pause; 

else
	
    load('Ch6_output.mat','kbarinitial', 'lbarinitial', 'averagehoursinitial', 'trbarinitial','taubinitial','tauwinitial');
    load('Ch6_output.mat','sp1initial', 'gninitial', 'massinitial', 'gbar2015', 'gainitial');
    load('Ch6_output.mat','loptinitial','coptinitial','aoptinitial','valueinitial');
    gbar = gbar2015;
end


% -------------------------------------------------------------------------
%
% STEP 4: Computation of the Final Steady State in year=year1
%
%--------------------------------------------------------------------------


calib=0;
if case_ss==1
    year1=year_final;				% 2010 or 2050: population parameters, multiple of 5
    year1=2050;
    if policy_change==1
        pen_repl=pen_repl_new;
    end

    clc
    ssfinal=1;
    disp('computation of the final steady state');
    gn=gnall(ntrans_finalss);
    disp('population growth rate final steady state: '); 
    gn
  %  pause;

    sp1=zeros(nage,1);
    sp1trans=sp1;
    for isp = 1:1:nage
        sp1(isp)=surviveprob(isp,ntrans_finalss);
        sp1trans(isp)=surviveprob(isp,ntrans);
    end

    speriod = linspace(20,20+nage-1,nage);
    figure
    plot(speriod,sp1initial,speriod,sp1,speriod,sp1trans);
    ylabel('Survival probabilities');
    xlabel('Age');
    legend('2015,2050,2100');
 %   pause;



    % computation of cohort mass in year_change 
    mass=ones(nage,1);
    for i=2:1:nage
        mass(i)=mass(i-1)*sp1(i-1)/(1+gn);
    end
    mass=mass/sum(mass);
    massfinal=mass;

    disp('dependency ratio final steady state: ');
    sum(mass(t+1:nage))/sum(mass(1:t))
    pause;

    trbar=trbarinitial;
    x0=[1.14532; 0.23909; 0.306724; 0.01788; 0.109246; 0.2480];
    disp('getss(x0)');
    y=getss(x0);
    y
  %  pause;
    xagg = fsolve(@(x)getss(x),x0);
    xagg = x0;
    disp('final steady state');
    xagg
    clc;
    disp('credit-constrained: ');
    fa(1) 
%    pause;

    xfinal=xagg;
    kbarfinal=xagg(1);
    lbarfinal=xagg(2);
    averagehoursfinal=xagg(3);
    trbarfinal=xagg(4);
    taubfinal=xagg(5);
    tauwfinal=xagg(6);
    sp1final=sp1;
    gnfinal=gn;

    gafinal=ga;
    loptfinal=lopt; 
    coptfinal=copt; 
    aoptfinal=aopt; 
    massfinal=mass;
    save('Ch6_output.mat','kbarfinal','lbarfinal','averagehoursfinal', 'trbarfinal', 'taubfinal', 'tauwfinal', '-append');
    save('Ch6_output.mat','sp1final', 'gnfinal', 'massfinal', 'gafinal', '-append');
    save('Ch6_output.mat','loptfinal', 'coptfinal', 'aoptfinal', '-append');


    ss_final=0;

else
    load('Ch6_output.mat','kbarfinal','lbarfinal','averagehoursfinal', 'trbarfinal', 'taubfinal', 'tauwfinal');
    load('Ch6_output.mat','sp1final', 'gnfinal', 'massfinal', 'gafinal');
    load('Ch6_output.mat','loptfinal', 'coptfinal', 'aoptfinal');
	
end

