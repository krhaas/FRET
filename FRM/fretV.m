function [V,F,Vpp,Fpp] = fretV(xref)

Vpp = pchip([0.5 0.8 0.9 1 1.1 1.2 1.5],[60 0 4 log(4) 5 log(2) 60]);
Fpp = fnder(Vpp);
Fpp.coefs = -Fpp.coefs;

V = ppval(Vpp,xref);
F = ppval(Fpp,xref);

return


tic
W1 = (4/7)*exp(-(xref-0.8).^2/2/0.001)/sqrt(2*pi*0.001);
W2 = (1/7)*exp(-(xref-1).^2/2/0.0005)/sqrt(2*pi*0.0005);
W3 = (2/7)*exp(-(xref-1.2).^2/2/0.00075)/sqrt(2*pi*0.00075);

p0 = W1+W2+W3;
V=-log(p0);
F=(-(xref-0.8).*W1/0.001-(xref-1).*W2/0.0005-(xref-1.2).*W3/0.00075)./p0;
toc

return


Vscale=15.47;
xscale=13.9;

x=x-1;
x=x*xscale;
V=(x.^6-15*x.^4+53*x.^2+2*x-15)/Vscale;
F=-(6*x.^5-60*x.^3+106*x+2)/Vscale*xscale;
dF=-(30*x.^4-180*x.^2+90)/Vscale*xscale*xscale;

return



Vscale=30;
xscale=15;

x=x-1;
x=x*xscale;
V=(x.^6-15*x.^4+50*x.^2+5*x-15)/Vscale;
F=-(6*x.^5-60*x.^3+100*x+5)/Vscale*xscale;
dF=-(30*x.^4-180*x.^2+90)/Vscale*xscale*xscale;
%V=V-min(V)-10;


%{
Vscale=30;
xscale=15;
x=x-1;
x=x*xscale;
V=(x.^6-15*x.^4+45*x.^2+5*x-15)/Vscale;
F=-(6*x.^5-60*x.^3+90*x+5)/Vscale*xscale;
dF=-(30*x.^4-180*x.^2+90)/Vscale*xscale*xscale;
%V=V-min(V)-10;
%}