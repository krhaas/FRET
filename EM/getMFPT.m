%Calculate Mean first passage time
function [kab,kba] = getMFPT(xA,xB,xref,V,D)

int1 = cumtrapz(xref,exp(-V));

Tab = (1/D)*integral(@(x)  exp(interp1(xref,V,x)).*...
    interp1(xref,int1,x) , xA, xB );

intall = trapz(xref,exp(-V));

Tba = (1/D)*integral(@(x)  exp(interp1(xref,V,x)).*...
    interp1(xref,intall-int1,x) , xA, xB );

kab = 1/Tab;
kba = 1/Tba;
end

