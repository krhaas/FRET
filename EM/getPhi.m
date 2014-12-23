function [py,gradD,dPhi,Phi,AB,LQ2out]=getPhi(Fdata,p0,CD,Chebs,pyNorm,LQ2)
%Takes in FRET data and model of dynamics and returns log-loss which is
%equal to negaitve of log-likelihood

x=Chebs{1}; w=Chebs{2}; Ldx=Chebs{3}; Tndx=Chebs{4}; Fndx=Chebs{12};
Vnx=Chebs{5}; Fnx=Chebs{6}; xref=Chebs{7}; w0=Chebs{8}; Np=Chebs{9}; Nh=Chebs{10};
Dnx=Chebs{11};

py=zeros(length(Fdata),1);
needgrad=nargout;

%D=exp(Dnx*CD);
D=Dnx*CD;

%iterate through each photon trajectory
    
Nv=64;%number of eignevectors of hamiltonian to compute and iterate.  increase for accuracy

%Tad=Fdata.Tad;%both acceptor and donor arrival times
iTad=Fdata.iTad;% =1 if acceptro and =0 if donor
Times=Fdata.Times;
DTs=Fdata.DTs;

IbA=Fdata.IbA;
IbD=Fdata.IbD;
Ba=Fdata.Ba;
Bd=Fdata.Bd;
IoA=IbA-Ba;
IoD=IbD-Bd;

Zofx=1./(1+xref.^6);
Iofx=(IoD*(1-Zofx)+Bd+IoA*Zofx+Ba);%dark intensity

N=length(xref);

p0=p0/sum(p0.*w0);
p02=p0;
p02=max(p0,1e-10);  %1e-10 Threshold so you only need 64 eigenvectors to get derivative right
%p02=max(p0,1e-20);
p02=p02/sum(p02.*w0);


Ke=spalloc(N,N,Np*Np*Nh);
for k=0:Nh-1
    range=1+k*(Np-1):k*(Np-1)+Np;
    Ke(range,range)=Ke(range,range)+...
        Ldx'*diag(D(range).*w.*p02(range))*Ldx;
end

Se=spdiags(w0.*p0,0,N,N);


%% Neumann No Ellimination of BCS
assert(issparse(Ke),'Ke isnt sparse')
warning('off','MATLAB:nearlySingularMatrix')
clear opts
opts.issym=1;
opts.v0=ones(size(Se,1),1);
opts.disp=0;
%opts.tol=eps/1e115;
try
    [Q,LQ,flag]=eigs(Ke,Se,Nv,'sm',opts);%generalized eigenvalue probelm, give smallest vec/vals
    assert(sum(flag) == 0,'nonconverged eigenvalues')
catch
    [Q,LQ,flag]=eigs(Ke,Se,Nv,-1,opts);%generalized eigenvalue probelm, give smallest vec/vals
    assert(sum(flag) == 0,'nonconverged eigenvalues')
end

    
LQ=diag(LQ);
%possible failures of numerical eigenvector method
assert(sum(imag(LQ)) == 0 , 'Imaginary Eigenvalues' ) 
LQ2out=LQ(2);

%smallest to largest rearrange
[LQ,iLQ]=sort(LQ);
Q=Q(:,iLQ);
LQ=LQ(1:Nv);
LQ(1)=0;
Qx=Q(:,1:Nv);
Qdx=Tndx*Qx;

%First Backward-Kolomogrov eignevector is flat
Qx(:,1)=1;  
Qdx(:,1)=0;
Qx2=Qx;     


%get actual eigenvector by multiplying by peq^(1/2)
Qx=spdiags(sqrt(p0),0,N,N)*Qx;

%Jeffery's prior
piA = 6*IoA*xref.^5./(1+xref.^6)./(IoA+Ba*(1+xref.^6));
piD = 6*IoD*xref.^5./(1+xref.^6)./(IoD*xref.^6+Bd*(1+xref.^6));
Prior = sqrt(piA.^2+piD.^2);

beta_a=(IoA+Ba)/Ba;
beta_d=(IoD+Bd)/Bd;

Ibeta_a=IoA/(1-1/beta_a);
Ibeta_d=IoD/(1-1/beta_d);

dt=1;

Jprior=36*xref.^10./(1+xref.^6).^3.*...
        (Ibeta_d*dt*(1-1/beta_d)^2./(xref.^6+1/beta_d)+Ibeta_a*dt*(1-1/beta_a)^2./(1+xref.^6/beta_a));

%Apply Dark to get diagnoalized H+Dark operator
Kdark = diag(LQ) + ...
    Qx'*spdiags(w0.*(Iofx-(Iofx)/LQ2.*log(Jprior)),0,N,N)*Qx;

[Qdark,LQdark]=eig(Kdark);
LQdark=diag(LQdark);
[LQdark,iLQdark]=sort(LQdark);
Qdark=Qdark(:,iLQdark);
Qxdark=Qx*Qdark;
Qx2dark=Qx2*Qdark;


%used in getting expectation of alpha and beta
iLQdark=zeros(Nv,Nv);
for i=1:Nv
    iLQdark(i,:)=1./(LQdark-LQdark(i));
    iLQdark(i,i)=0;
end

%get acceptor/donor photon operators
PA=Qxdark'*spdiags(w0.*(IoA*Zofx+Ba),0,N,N)*Qxdark;
PD=Qxdark'*spdiags(w0.*(IoD*(1-Zofx)+Bd),0,N,N)*Qxdark;


%initialize with equilibrium distribution
btemp=zeros(Nv,1);
btemp(1)=1;
btemp=Qdark'*btemp;

%Prior probability Initialization 

%piAD = sqrt(piA.*piD);
%btemp = Qxdark'*piAD;
atemp=btemp;
%Apply initial wait time to start and end of trajectory


%% ALPHA RECURSION 
%Run alpha Mex code that iterates through photons for alpha
[Anorm,alpha] = runAlpha(atemp,DTs,iTad,PA,PD,LQdark);

Anorm(Times+2)=sum(btemp.*atemp);%make sure integral norms to 1
atemp=atemp/Anorm(end);
btemp=btemp/Anorm(Times+2)/Anorm(Times+1);
%betax(:,end)=Qx*btemp;

%% Calculate Likelihood 
py=-sum(log(Anorm)+pyNorm/(Times+2));

%if need gradient keep going
if ( needgrad > 1 )
    
%% BETA Recursion    
%Flip iLQdark because MATLAB outputs column major and C in row major
G = runBeta(btemp,DTs,iTad,PA,PD,LQdark,alpha,Anorm,iLQdark');
G = G';

%Calculate Initial and final alpha & beta
AB=(Qx2dark*(atemp+btemp));  

%Transfer from Dark basis to pure H basis
GHo=Qdark*G*Qdark'; 

%Calculate relative diffusion derivative
Pdxdx = sum(Qdx.*(Qdx*GHo),2);

gradD=sum(w0.*p0.*Pdxdx);

% Include below for infinite time start/end conditions
% Ginit will just be zero everywhere except Ginit(1,2:end)
Ginit=Qdark*(btemp+atemp)./LQ;
Ginit(1)=0;
Ginit=[1 ;zeros(Nv-1,1)]*Ginit';
%GHo = GHo + Ginit;
%}

%Calculate relative time averaged position 
Phi = sum(Qx2.*(Qx2*GHo),2);

%Calculate relative force Derivative
dPhi = sum(Qx2.*(Qdx*GHo),2) + sum(Qdx.*(Qx2*GHo),2);

gradF=-D.*p0.*dPhi/2;

end
end


%{
% Old Alpha
Anorm=ones(1,Times+2);
Anorm(1)=norm(atemp);
%atemp=atemp/Anorm(1);
alpha=zeros(Nv,Times+1);
alpha(:,1)=atemp;

%alphax(:,1)=Qx*atemp;

%iterate through trajectory for alphas
for i=1:Times
    %DTs=Tad(i)-Tad(i-1);%time between photons
    
    if( iTad(i) ) %apply donor or acceptor operator
        atemp = PD*(atemp.*exp(-LQdark*DTs(i)));
    else
        atemp = PA*(atemp.*exp(-LQdark*DTs(i)));       
    end
    
    Anorm(i+1)=norm(atemp);
    atemp=atemp/Anorm(i+1);
    alpha(:,i+1)=atemp;
   
    %alphax(:,i)=Qx*atemp;
end
%}

%{ 
    %Old Beta Recursion
    G=zeros(Nv,Nv);%gradient expectation
someones=ones(1,Nv);
for i=Times:-1:1
    %DT=Tad(i+1)-Tad(i);
    DT=DTs(i);
    expLQDT=exp(-LQdark*DT);
    
    if( iTad(i) )
        btemp = PD*btemp;
    else
        btemp = PA*btemp;
    end
    
    %time smooth the derivative of operator
    eLQDT=(expLQDT*someones);
    outer=(alpha(:,i)*btemp');
    G=G+outer.*(diag(expLQDT*DT)+iLQdark.*(eLQDT-eLQDT'));
    
    btemp=btemp.*expLQDT;
   
    %beta(:,i)=btemp;
    %betax(:,i)=Qx*btemp;
    btemp=btemp/Anorm(i);

end
%}


%{ 
Old version of PxDx stuff

Pdxdx=zeros(size(xref));
for i=1:Nv
    for j=1:Nv
        Pdxdx=Pdxdx+Qdx(:,i).*Qdx(:,j).*GHo(i,j);
    end
end

%Calculate statistic of force Derivative using modified vector
Pxdxdxx=zeros(size(xref));
for i=1:Nv
    for j=1:Nv
        Pxdxdxx=Pxdxdxx+(Qx2(:,i).*Qdx(:,j)+Qx2(:,j).*Qdx(:,i)).*GHo(i,j);      
    end
end

Pxx=zeros(size(xref));
for i=1:Nv
    for j=1:Nv
        Pxx=Pxx+Qx2(:,i).*Qx2(:,j).*GHo(i,j);
    end
end

%{
for i=1:Nv
    for j=1:Nv
        Pxdxdxx=Pxdxdxx+(Qx2(:,i).*Qdx(:,j)+Qx2(:,j).*Qdx(:,i)).*Ginit(i,j);
        Pxx=Pxx+Qx2(:,i).*Qx2(:,j).*Ginit(i,j);
    end
end
%}


%}

%% Dirichlet Ellimination of BCS
%{
%BC1=Tndx(1,:);
%BC2=Tndx(end,:);
%BC=[BC1;BC2];
BC=zeros(2,N);
BC(1,1)=1;
BC(end,end)=1;
NullM=null(BC);
NullM=sparse(NullM);
Ke=NullM'*sparse(K)*NullM;
Se=NullM'*Se*NullM;
warning('off','MATLAB:nearlySingularMatrix')
clear opts
opts.issym=1;
opts.v0=ones(size(Se,1),1);
%opts.disp=0;
opts.tol=eps/1e50;
[Qe,LQ]=eigs(Ke,Se,Nv,'sm',opts);
Q=NullM*Qe;
%}
%%