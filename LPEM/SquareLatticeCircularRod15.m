function [xt,epst]=SquareLatticeCircularRod15(d,r,epsrod,epssub,n)

% This function gives profile of a unit cell of square lattice with
% circular rods:
%                        2*r
%                   <--------->
%        -     
%       |      
%       |               ----
%       |            /********\
%       |           |**********|
%    d <            |**********|
%       |            \********/
%       |               ----
%       |      
%        -                   
%              <------------------->
%               Grating period (=1)
% 
% d : Grating thickness (normalized to grating period)
% r : Rod's radius (normalized to grating period)
% epsrod : Permitttivity of rods
% epssub : Permitttivity of substrate
% n : L*nlayer


z=LPEM_zgen15(d,n);
zlen = length(z);

xt=zeros(2,zlen);
epst = zeros(3,zlen);

epst(1,:)= epssub;
epst(2,:) = epsrod;
epst(3,:) = epssub;

xt(1,:)=.5;
xt(2,:)=.5;

for k=1:zlen
    if (z(k)>d/2-r)&&(z(k)<d/2+r)
        xt(1,k)=1/2-sqrt(r^2-(z(k)-d/2)^2);
        xt(2,k)=1/2+sqrt(r^2-(z(k)-d/2)^2);
    end
end
