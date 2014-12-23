function [x,w,dL,L]=getLagrange(N,xratio,varargin)

N=N-1;
[x,w,P]=lglnodes(N);

x=flipud(x);
w=flipud(w);
Pend=flipud(P(:,end));


dL=zeros(size(P));
for i=1:N+1
    dL(:,i)=Pend./Pend(i)./(x-x(i));
    dL(i,i)=0;
end
dL(1,1)=-(N*(N+1))/4;
dL(end,end)=N*(N+1)/4;
L=diag(x);
if( length(varargin) ~= 0 )
 
    xref = varargin{1};
    xref=xref-mean(xref);
    xref=xref/xratio;
    L=ones(length(xref),N+1);
    dL=zeros(length(xref),N+1);
    for i=1:N+1
        for j=1:N+1
            if( i ~= j)
                L(:,i)=L(:,i).*(xref-x(j))/(x(i)-x(j));
            end
        end
    end
        
   
    for i=1:N+1
        
        for k=1:N+1;
            if( k ~= i )
           temp=xref*0+1;
        for j=1:N+1
            if( i ~= j && j ~=k )
                temp=temp.*(xref-x(j))/(x(i)-x(j));
            end
            if( i ~= j && j == k )
                temp=temp/(x(i)-x(j));
            end
        end
        
        dL(:,i)=dL(:,i)+temp;
            end
        end
    end
    
    
end

x=x*xratio;
w=w*xratio;
dL=dL/xratio;

return
P=zeros(length(xref),N+1);
    P(:,1)=1;    P(:,2)=xref;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*xref.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
    dPend=(N+1)*(xref.*P(:,N+1)-P(:,N))./(xref.^2-1);
    dPend(end)=N*(N+1)/2;
    dPend(1)=(-1)^(N+1)*N*(N+1)/2;
    PendBig=P(:,end);
        
    dL=zeros(size(P));
    for i=1:N+1
       dL(:,i)=PendBig./Pend(i)./(xref-x(i));
       dL(i,i)=0;
    end
    dL(1,1)=-(N*(N+1))/4;
    dL(end,end)=N*(N+1)/4;

    L=zeros(size(P));
    for i=1:N+1;
        L(:,i)=(xref.^2-1)./N./(N+1)./Pend(i)./(xref-x(i)).*dPend;
    end