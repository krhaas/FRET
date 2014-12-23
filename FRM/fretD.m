function [D,dD] = fretD(x)
    D=(1+0.5*tanh((x-1)*10))*4e2;
    dD=(0.5*10*(sech((x-1)*10)).^2)*4e2;
end
