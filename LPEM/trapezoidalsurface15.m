function [xt,epst,xt_original,z]=trapezoidalsurface15(d1,d2,d3,d4,eps1,eps3,d,n)
% This function gives profile of a trapezoidal surface relief grating.
% d1, d2, d3 and d4 : Place of angles of trapezoid in one period and on the
% x axis (0<=di<=1 , i=1,2,3,4)
%                  ---------
%                /           \
%               /             \
%              /               \
%              -----------------
%             d1  d2       d3  d4
% eps1 : Permittivity of incident region
% eps3 : Permittivity of lower region
% d : grating thickness (normalized to grating period)
% n : L*nlayer
z=LPEM_zgen15(d,n);
zlen = length(z);

xt=zeros(2,zlen);
epst = zeros(3,zlen);

epst(1,:)= eps1;
epst(2,:) = eps3;
epst(3,:) = eps1;

for k=1:zlen
    xt(1,k)=(d1-d2)*z(k)/d+d2;
    xt(2,k)=(d4-d3)*z(k)/d+d3;
    
    xt_original(1,k)=xt(1,k);
    xt_original(2,k)=xt(2,k);
    
    %%{
    if xt(1,k)<1 && xt(2,k)>1
        x0=xt(1,k);
        xt(1,k)=mod(xt(2,k),1);
        xt(2,k)=x0;
        epst(1,k)= eps3;
        epst(2,k) = eps1;
        epst(3,k) = eps3;
    elseif xt(1,k)>1
        xt(1,k)=mod(xt(1,k),1);
        xt(2,k)=mod(xt(2,k),1);
    end
    %}
        
end
