function [xt,epst]=TriangularLatticeCircularRodDefect15(d,r,rdefect,epsrod,epssub,epsroddefect,epssubdefect,n)
 
% This function gives profile of a unit cell of triangular lattice with
% circular rods and the central rod's radius is different:
%                r
%              <---->
% 
%        -     *****|         |*****
%       |      ****/           \****
%       |      --                 --
%       |            2*rdefect
%       |              <---->          -
%       |                --             |
%       |              /****\            > defect
%    d <               \****/           | 
%       |                --             |
%       |                              -
%       |
%       |      --                 --
%       |      ****\           /****
%        -      *****|         |*****  
%              
%              <------------------->
%               Grating period (=1)
% 
% d : Grating thickness (normalized to grating period)
% r : Rod's radius (normalized to grating period)
% rdefect : Radius of defect's rod (normalized to grating period)
% epsrod : Permittivity of defect's rod
% epssub : Permittivity of defect's substrate
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
    epst(:,k)=[epssub epssub epssub epssub epssub];
    if z(k)<r
        xt(1,k)=sqrt(r^2-z(k)^2);
        xt(4,k)=1-xt(1,k);
        epst(:,k)=[epsrod epssub epsrod epssub epsrod];
    end
    if (z(k)>d/2-rdefect)&&(z(k)<d/2+rdefect)
        xt(2,k)=1/2-sqrt(rdefect^2-(z(k)-d/2)^2);
        xt(3,k)=1/2+sqrt(rdefect^2-(z(k)-d/2)^2);
        epst(:,k)=[epsroddefect epssubdefect epsroddefect epssubdefect epsroddefect];
    end
 
    if z(k)>d-r
        xt(1,k)=sqrt(r^2-(z(k)-d)^2);
        xt(4,k)=1-xt(1,k);
        epst(:,k)=[epsrod epssub epsrod epssub epsrod];
    end
end