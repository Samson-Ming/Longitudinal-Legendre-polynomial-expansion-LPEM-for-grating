function [xt,epst,xt_original,z]=overhang_sinusoidalsurface15(eps1,eps3,d,Phi,n)

% This function gives profile of a sinusoidal surface relief grating.
% 
%     --          --          --
%   /    \      /    \ eps1 /
%  /      \    / eps3 \    /
%-          --          -- 
%            <---------->
%          Grating period (=1)
% 
% eps1 : Permittivity of incident region
% eps3 : Permittivity of lower region
% d : grating thickness (normalized to grating period)
% n : L*nlayer

z=LPEM_zgen15(d,n);
zlen = length(z);

xt=zeros(2,zlen);
xt_original=zeros(2,zlen);
epst = zeros(3,zlen);

epst(1,:)= eps1;
epst(2,:) = eps3;
epst(3,:) = eps1;


for k=1:zlen
    xt(1,k)=acos(-(2*(d-z(k))/d-1))/(2*pi/cos(Phi));
    xt(2,k)=cos(Phi)-xt(1,k);
    
    xt(1,k)=sec(Phi)*xt(1,k)+tan(Phi)*(d-z(k));
    xt(2,k)=sec(Phi)*xt(2,k)+tan(Phi)*(d-z(k));
    
    xt_original(1,k)=xt(1,k);
    xt_original(2,k)=xt(2,k);
    
    %{
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
    
    %%{
    x1=mod(xt(1,k),1);
    x2=mod(xt(2,k),1);
    if x1>x2
        xt(1,k)=x2;
        xt(2,k)=x1;
        epst(1,k)= eps3;
        epst(2,k) = eps1;
        epst(3,k) = eps3;
    else
        xt(1,k)=x1;
        xt(2,k)=x2;
    end
%}
    
end
