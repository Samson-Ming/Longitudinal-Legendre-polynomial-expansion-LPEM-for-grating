function [xt,epst]=TriangularLatticeCircularRod15(d,r,epsrod,epssub,n)

% This function gives profile of a unit cell of triangular lattice with
% circular rods:
%                r
%              <---->
% 
%        -     *****|         |*****
%       |      ****/           \****
%       |      --                 --
%       |               ----
%       |            /********\
%       |           |**********|
%    d <            |**********|
%       |            \********/
%       |               ----
%       |      --                 --
%       |      ****\           /****
%        -     *****|         |*****  
%              
%              <------------------->
%               Grating period (=1)
% 
% d : Grating thickness (normalized to grating period)
% r : Rod's radius (normalized to grating period)
% epsrod : Permittivity of rods
% epssub : Permittivity of substrate
% n : L*nlayer

z=LPEM_zgen15(d,n);
zlen = length(z);

xt=zeros(4,zlen);
epst=zeros(5,zlen);

epst(1,:)=epsrod;
epst(2,:)=epssub;
epst(3,:)=epsrod;
epst(4,:)=epssub;
epst(5,:)=epsrod;

xt(1,:)=0;
xt(2,:)=.5;
xt(3,:)=.5;
xt(4,:)=1;
for k=1:zlen
    if z(k)<r
        xt(1,k)=sqrt(r^2-z(k)^2);
        xt(4,k)=1-xt(1,k);
    end
    if abs(z(k)-d/2)<r
        xt(2,k)=1/2-sqrt(r^2-(z(k)-d/2)^2);
        xt(3,k)=1-xt(2,k);
    end
    if abs(z(k)-d)<r
        xt(1,k)=sqrt(r^2-(z(k)-d)^2);
        xt(4,k)=1-xt(1,k);
    end
end
