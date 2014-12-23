%function runPhiEigEM(varargin)
% Get EM optimized smFRET 
% runPhiEigEM(smFRETdata,etaF,run#)
%% Run as Script
FRETname='QuadFdatax2.mat';
etaF=2e-7;
run=3;
loops=10001;
recordIter = 50;

if( nargin < 4 )
    'Need to pass name of FRET trajectory and etaF'
else
    FRETname = varargin{1}
    etaF = str2num(varargin{2})
    run = str2num(varargin{3})
    loops = str2num(varargin{4})
end

format long

load(FRETname)


%select trajectory number to work on.
Fdata=Fdata(run);
Fdata.Times=50000; %Number of photons total
Tad=Fdata.Tad(1:Fdata.Times);
Fdata.iTad=(Fdata.iTad(1:Fdata.Times) > eps );
Fdata.DTs = Tad - [0; Tad(1:end-1)];


%[xref,Tnx,Tndx,w0,Tnix,Tnddx]=getCheby(Res,N,xratio);

Nh=256;%Spectral element number
Np=8;%Spectral element order
N=Nh*Np-Nh+1;

meanPos=sum(Fdata.xc(:,1).*Fdata.xc(:,2))/sum(Fdata.xc(:,2));
stdevPos=sqrt(sum((Fdata.xc(:,1)-meanPos).^2.*Fdata.xc(:,2))./sum(Fdata.xc(:,2)));

maxx=max(Fdata.xc(:,1));
minx=min(Fdata.xc(:,1));

maxx=Fdata.xref(find(Fdata.rawpdf > 1e-2,1,'last'));
minx=Fdata.xref(find(Fdata.rawpdf > 1e-2,1,'first'));


xratio=(maxx-minx)/2;% Lenght L of box
xavg=(maxx+minx)/2;% Center of Box x=r/Ro
xratio=0.6;
xavg=1.0;


[x,w,Ldx]=getLagrange(Np,xratio/Nh);

xref=zeros(N,1);% grip points of domain
w0=zeros(N,1);% integration kernel (int=sum(w0.*f(x)))
Tndx=spalloc(N,N,Np*Np*Nh);%derivative of basis functions
for i=0:Nh-1
    xref(1+i*(Np-1):i*(Np-1)+Np)=x+2*xratio/Nh*i+xratio/Nh;
    w0(1+i*(Np-1):i*(Np-1)+Np)=w0(1+i*(Np-1):i*(Np-1)+Np)+w;
    Tndx(1+i*(Np-1):i*(Np-1)+Np,1+i*(Np-1):i*(Np-1)+Np)=Tndx(1+i*(Np-1):i*(Np-1)+Np,1+i*(Np-1):i*(Np-1)+Np)+Ldx;
end
for i=0:Nh-1
    Tndx(Np+i*(Np-1),:)=Tndx(Np+i*(Np-1),:)/2;
end
Tndx(end,:)=Tndx(end,:)*2;

xref=xref+xavg-xratio;


Dnx=xref*0+1;

L=(xref(end)-xref(1));
p0=cos((xref-xavg)*pi/L);
p0(1)=0; p0(end)=0;
dp0=-sin((xref-xavg)*pi/L)*pi/L;

pref = exp(-fretV(xref));
pref = pref/sum(w0.*pref);

Chebs{1}=x; Chebs{2}=w; Chebs{3}=Ldx; Chebs{4}=Tndx; Chebs{11}=Dnx; Chebs{12}=0; Chebs{14}=dp0;
Chebs{5}=0; Chebs{6}=0; Chebs{7}=xref; Chebs{8}=w0; Chebs{9}=Np; Chebs{10}=Nh; Chebs{13}=p0; 


restart=1;
if(restart == 1 )
    
pCP=p0; %Cutting planes
dpCP=dp0;
%pCP=p0.*xref;
%dpCP=p0+dp0.*xref;

dpCP=dpCP/sqrt(sum(w0.*pCP.^2));
pCP=pCP/sqrt(sum(w0.*pCP.^2));


pB=pCP; %Optimal point
dpB=dpCP;
etaD=NaN;
pyCPs=zeros(1,ceil(loops/recordIter));
%epsilon=pyCPs;
pyBs=pyCPs;
%FreeEs=pyCPs;
%ldCPs=pyCPs;
ldBs=pyCPs;
etaDs=pyCPs;
etaFs=pyCPs;
kabs=pyCPs;
kbas=pyCPs;
Iteration=pyCPs;
%alphas=zeros(loops,loops);
%gradDs=pyCPs;
%pCPs=zeros(length(xref),loops);
%dpCPs=pCPs;
%ABs=pCPs;
%Phis=pCPs;
%dPhis=pCPs;

pBs=zeros(length(xref),ceil(loops/recordIter));
dpBs=pBs;
%epsilon(1)=NaN;

%[pyB,gradDB,dPhi,Phi] = getPYDPYBundleEM(Fdata,pref,500,Chebs,0);
%[pyB2,gradDB2,dPhi2,Phi2] = getPYDPYBundleEMold(Fdata,pref,500,Chebs,0);

%% Optimize initial Diffusion value
%pyNorm=getPYDPYBundleEM(Fdata,pB.^2,50,Chebs,0);
tauI=0;
options4 = optimset('Display','iter-detailed','GradObj','on','TolFun',1e-4,...
        'LargeScale','off','DerivativeCheck','off','FinDiffType','central','FinDiffRelStep',1e-3);
[ldB,pyNorm,exitflag,output,grad,hessian] = fminunc(@(ldB) getPhi(Fdata,pB.^2,ldB,Chebs,0),50,options4);

%pyNorm=getPYDPYBundle2(Fdata,pB.^2,ldB,Chebs,0);
else
    %pCP=pCPs(:,restart);
    %dpCP=dpCPs(:,restart);
    pB=pBs(:,restart);
    dpB=dpBs(:,restart);
    etaD=etaDs(restart);
    etaF=etaFs(restart);
    %alpha=alphas(1:restart,restart);
    ldCP=ldCPs(restart);
    ldB=ldBs(restart);
end

record=1;
for iter=restart:loops
    iter
    pB = pCP;
    dpB = dpCP;
    
    %Goto Optimum Diffusion point
    if( ~mod(iter,5) )
    options4 = optimset('Display','iter-detailed','GradObj','on','TolFun',1e-4,...
        'LargeScale','off','InitialHessType','user-supplied','InitialHessMatrix',hessian);
    [ldB,pyB,exitflag,output,grad,hessian] = ...
        fminunc(@(ldB) getPhi(Fdata,pB.^2,ldB,Chebs,pyNorm),ldB,options4);
    etaD = ldB*etaF*sum(w0.*dpB.^2);
    end
    
    
    [pyB,gradDB,dPhi,Phi,AB]=getPhi(Fdata,pB.^2,ldB,Chebs,pyNorm);

   
    pstar = (pB.^2.*Phi).^(1/(1+ldB*etaF));
    pstar = pstar/sum(w0.*pstar);
    pCP = sqrt(pstar);
    dpCP = Tndx*pCP;
    
   %% Plot Iterations
    if( ~mod(iter-1,recordIter) )
        
    [kab,kba] = getMFPT(0.8,1.2,xref,-log(pB.^2),ldB);
    kabs(record)=kab;
    kbas(record)=kba;
    pyBs(record)=pyB; 
    etaDs(record)=etaD;
    etaFs(record)=etaF;
    
    %FreeEs(record)=pyB + etaF*ldB*sum(w0.*dpB.^2);
    pBs(:,record)=pB;
    dpBs(:,record)=dpB;
    ldBs(record)=ldB;
    Iteration(record)=iter;
    
    clf
    set(gcf, 'PaperUnits','inches','PaperPosition', [0.5 0.5 7.5 10]);
    subplot(3,1,1)
    
    [AX,H1,H2] = plotyy(xref,[pB.^2 pCP.^2]...
        ,xref,[dPhi]);
    title(sprintf('Diffusion = %0.4G', ldB))
    set(AX(1),'XLim',[xref(1) xref(end)]);
    set(AX(1),'YLim',[0 10]);
    set(AX(1),'YTick',0:2:10);
    set(AX(2),'XLim',[xref(1) xref(end)]);
    ylabel('Probability')
    xlabel('Position (R/R_0)')
    set(get(AX(2),'Ylabel'),'string','\Phi\prime(x) Statistic')
    
    subplot(3,1,2)
    plot(Iteration(1:record),-pyBs(1:record),'-o');
    title(['{\itk}_{A\rightarrowB} = ' sprintf('%0.4G', kab) ...
        '   {\itk}_{B\rightarrowA} = '  sprintf('%0.4G',kba)])
    
    ylabel('Log-Likelihood')
    xlabel('Iteration #')
   	
    
    %legend('Log-Likelihood','\epsilon = J^*-J_C_P','Location','North','Orientation','horizontal')
    
    
    subplot(3,1,3)
    [AX] = plot(Iteration(1:record), (ldBs(1:record)),'-o');
    xlabel('Iteration #')
    title(['\eta_F' sprintf(' = %0.4G',etaF) '     \eta_D' sprintf(' = %0.4G',etaD)])
    ylabel('Diffusion Constant')
    
    drawnow
    print(gcf,'-dpdf',[pwd '/EMfig/' FRETname sprintf('etaF%.2dEM%.2d.pdf',etaF,iter)]);
    
    
    save([pwd '/EMdata/' FRETname sprintf('etaF%.2dEM%.2d.mat',etaF,iter)],...
        'ldB','pB','dpB','Phi','dPhi','xref','etaD','kab','kba') 
    
    
    record=record+1;
    end
   
   
end

save([FRETname sprintf('etaF%.2dEM.mat',etaF)],...
    'ldBs','pyBs','pBs','dpBs','Phi','dPhi','etaDs','etaFs','xref','kabs','kbas','Iteration')

%end
