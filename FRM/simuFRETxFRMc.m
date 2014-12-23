% Simulation of smFRET trajectory.

Nruns=12;
parfor runs=1:Nruns
dt=0.000000001;
kbT=1;
D=500;
totPhotons=100000;
alpha=0.1;
%dA=linspace(0.7,1.3,Nruns);

Ta=zeros(totPhotons,1);
Td=zeros(totPhotons,1);

aRecord=0;
dRecord=0;

xref=linspace(0,2,500);
%V1=fretV(xref);
%V1=2000*(1-xref).^4-200*(1-xref).^2+2*(1-xref);
%pxref=exp(-V1/kbT)/sum(exp(-V1/kbT));
[V,F,Vpp,Fpp] = fretV(xref);

cumPxref=cumtrapz(xref,exp(-V));
cumPxref=cumPxref/cumPxref(end);
xpos=xref(find(cumPxref>rand,1));


%Time by 10
%parameters
%res=100;

Id0 = 8000*res;
Ia0 = 15000*res;
Bd = 10*res;
Ba = 20*res;

%{
Old Parameters
Id0=4900*res;
Ia0=4900*res;
Bd=50*res;
Ba=50*res;
%}


beta_a=(Ia0+Ba)/Ba;
beta_d=(Id0+Bd)/Bd;

Ibeta_a=Ia0/(1-1/beta_a);
Ibeta_d=Id0/(1-1/beta_d);


t=0;
Idt=Id0*(xpos.^6./(1+xpos.^6))+Bd;
Iat=Ia0*(1./(1+xpos.^6))+Ba;

%random reaction likelihoods
FRETranda=-log(rand);
FRETrandd=-log(rand);
FRETcuma=0;
FRETcumd=0;

iter=1;
xsave=zeros(totPhotons,2);
while( sum(aRecord) + sum(dRecord) < totPhotons )
    
    %[V,F]=fretV(xpos);
    %[D,dD]=fretD(xpos);
    %D=500; dD=0;
    %F=F+dD/D;
    
    F = ppval(Fpp,xpos);
    
    zeta=(1/(1+xpos^6));
    Idt=Id0*(1-zeta)+Bd;
    Iat=Ia0*(zeta)+Ba;
    
    tstep=dt;
    %time they would hit
    tHITa=(FRETranda-FRETcuma)/Iat;
    tHITd=(FRETrandd-FRETcumd)/Idt;
    
    %do they hit in this time window
    
        
    if(tHITa < dt && tHITa<tHITd)
        FRETcuma=-Iat*tHITa;
        FRETranda=-log(rand);
        aRecord=aRecord + 1;
        Ta(aRecord)=t+tHITa;
        
        xsave(aRecord+dRecord,:)=[t+tHITa,xpos];
        tstep=tHITa;
        
    end
    
    if(tHITd < dt && tHITa>tHITd )
        FRETcumd=-Idt*tHITd;
        FRETrandd=-log(rand);
        dRecord=dRecord + 1;
        Td(dRecord)=t+tHITd;
        
        xsave(aRecord+dRecord,:)=[t+tHITd,xpos];
        tstep=tHITd;
    
    end
    
    FRETcuma=(FRETcuma+Iat*tstep);
    FRETcumd=(FRETcumd+Idt*tstep);
    
    
    t=t+tstep;
    xpos=xpos+D*tstep*F+sqrt(tstep*2*kbT*D)*randn;
    iter=iter+1;
    %xsave(iter,:)=[t,xpos];
end

%iter
%xsave=xsave(1:iter-1,:);
Fdata=struct;
Fdata.Xsave=xsave;

maxAi=find(Ta>0 ,1,'last');
maxDi=find(Td>0 ,1,'last');


Ta=Ta(1:maxAi);
Td=Td(1:maxDi);
Fdata.Ta=Ta;%Acceptor donor arrival times
Fdata.Td=Td;


Fdata.Ableach=Ta(maxAi);
Fdata.Dbleach=Td(maxDi);
Fdata.IbA=Ibeta_a;%intensity with background
Fdata.IbD=Ibeta_d;
Fdata.Ba=Ba;%background
Fdata.Bd=Bd;


[Tad,iTad]=sort([Ta;Td]);
maxA=length(Ta);
maxD=length(Td);
iTad=iTad>maxA;
iTad=ones(size(iTad)).*iTad;
Times=maxA+maxD;
Fdata.Tad=Tad;
Fdata.iTad=iTad;
Fdata.Times=Times;
Fdata.Alpha=alpha;
xc=zeros(length(Tad),4);
ixc=1;
start=0;
tb=1;
rawpdf=xref*0;
for tf=1:Times
    nd=sum(iTad(tb:tf) == 1);
    na=sum(iTad(tb:tf) == 0);
    dt=Tad(tf)-start;
    
    if( na == 0 || nd == 0)
        continue
    end

    xtest=(beta_a/beta_d*(Ibeta_d*na-Ibeta_a*nd*beta_d)/(Ibeta_a*nd-Ibeta_d*na*beta_a))^(1/6);
    Jtest=36*xtest^10/(1+xtest^6)^3*...
        (Ibeta_d*dt*(1-1/beta_d)^2/(xtest^6+1/beta_d)+Ibeta_a*dt*(1-1/beta_a)^2/(1+xtest^6/beta_a));
    atest=sqrt(1/Jtest);
    if( atest < alpha && imag(atest) == 0 )
        rawpdf=rawpdf+dt*exp(-(xtest-xref).^2/2/atest^2)/sqrt(2*pi*atest^2);
        xc(ixc,:)=[xtest,dt,Tad(tf)-dt/2,atest];
        ixc=ixc+1;
        tb=tf+1;
        start=Tad(tf);
    end
end
xc=xc(1:ixc-1,:);
Fdata.xc=xc;
rawpdf=rawpdf/trapz(xref,rawpdf);
Fdata.rawpdf=rawpdf;
Fdata.xref=xref;
Fdata.avgAlpha=sum(xc(:,4).*xc(:,2))/sum(xc(:,2));
%Fdata.D=fretD(xref);
Fdata.D = D;
Fdata.V=fretV(xref);

Fbig(runs)=Fdata;
end

Fdata=Fbig;

save(['QuadFdatax' sprintf('%d',res) '.mat'],'Fdata')
