%Calculate Mean first passage time
function [kab,kba] = getMFPT(xA,xB,xref,V,D)
%{
Nh=256;%Spectral element number
Np=8;%Spectral element order
N=Nh*Np-Nh+1;

xratio=0.5;
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


V = fretV(xref);
D = 500;
%}


int1 = cumtrapz(xref,exp(-V));

Tab = (1/D)*integral(@(x)  exp(interp1(xref,V,x)).*...
    interp1(xref,int1,x) , xA, xB );

intall = trapz(xref,exp(-V));


Tba = (1/D)*integral(@(x)  exp(interp1(xref,V,x)).*...
    interp1(xref,intall-int1,x) , xA, xB );

kab = 1/Tab;
kba = 1/Tba;
%ratio = kab/kba
end

