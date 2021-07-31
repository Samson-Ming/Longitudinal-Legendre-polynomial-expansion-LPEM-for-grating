function [xt,epst]=sinusoidalsurface15(eps1,eps3,d,n)

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
epst = zeros(3,zlen);

epst(1,:)= eps1;
epst(2,:) = eps3;
epst(3,:) = eps1;


for k=1:zlen
    xt(1,k)=acos(2*z(k)/d-1)/(2*pi);
    xt(2,k)=1-xt(1,k);
end
