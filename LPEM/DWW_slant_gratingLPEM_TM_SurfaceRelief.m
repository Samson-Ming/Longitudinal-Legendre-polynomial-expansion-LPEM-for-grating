% In The Name of GOD
% 
% This code has been written for analysis of longitudinally nonhomogenous
% volume gratings. The incident wave is TM polarized.
% 
% In the first section (Inputs) you should determine the structure and
% adjust convergence parameters. Here, an example of defining a structure
% is given.
% 
% Example: A trapezoidal surface relief grating 
% 
% [xt,epst]=trapezoidalsurface15(0,0.2,0.6,0.8,eps1,eps3,d,L*nlayer);
% 
% 
% Incident plane wave(TM)
%                   \    |
%                    \   |
%                     \  |   Region 1 (eps1)
%                      \ |
%                      ---         ---          -
%                     /   \       /   \          |
%                    /     \     /     \          > d
%                   /       \   /       \        |
%                ---         ---         ---    -
%                   <---------->            
%                 Grating Period(=1)         
%                  
%                            Region 3 (eps3)
% 
% The outputs of this program are diffraction efficiencies which are stored
% in vectors 'DE1' and 'DE3' for reflected and transmitted waves,
% respectively. For example for N = 2 we have:
% DE1 = diffraction efficiency of [2 1 0 -1 -2]'th reflected order
% DE3 = diffraction efficiency of [2 1 0 -1 -2]'th transmitted order
 
%%%%%%%
%INPUTS
%%%%%%%
clear
clc
% ----------------------
% Convergence parameters
% ----------------------
% It is recommended not to change "M" and "nlayer".
N=30;     %Truncation order
L=40;       % Number of layers in a unit cell.
M=3;        % Number of polynomial terms in each layer.
nlayer = 1;  % this parameter increase the integration accuracy
% ---------------------------
% Global structure parameters
% ---------------------------
%grating period must be normalized to 1
lambda=0.525/0.425;     %Vacuum wavelength (Normalized to the grating period)
d=0.36/0.425;      % Grating thickness (Normalized to the grating period)
eps1=1;     %Permittivity of incident medium
%eps3=1.72^2;     %Permittivity of transmission medium
eps3=(1-5i)^2;     %Permittivity of transmission medium
alpha=15*pi/180;% Angle of incidence (in radian)
%alpha=0*pi/180;% Angle of incidence (in radian)
% ------------------------
% Specifying surface type
% ------------------------
% Two different surfaces can be defined here: sinusoidal and trapezoidal
% surfaces.
%Phi=0*(pi/180);
Phi=15*(pi/180);
ff=0.5;  %filling factor
[xt,epst]=trapezoidalsurface15(0,d*tan(Phi),d*tan(Phi)+ff,ff,eps1,eps3,d,L*nlayer);
% [xt,epst]=sinusoidalsurface15(eps1,eps3,d,L*nlayer);
 
% -------------------
% initial computations
% -------------------
k0=2*pi/lambda;
k1=k0*(eps1^.5);
k3=k0*(eps3^.5);
ux=cos(alpha);
uz=-sin(alpha);
kx = k1 * sin(alpha) - (-N:N) * 2 * pi;
kz1 = sqrt(k1^2-kx.^2);
kz3 = sqrt(k3^2-kx.^2);
for j=1:2*N+1
    if real(kz1(j))>0
       kz1(j)=-kz1(j);
    end
    if imag(kz1(j))<0
       kz1(j)=-kz1(j);
    end
    if real(kz3(j))<0
       kz3(j)=-kz3(j);
    end
    if imag(kz3(j))>0
       kz3(j)=-kz3(j);
    end
end
Kz1=diag(kz1);
Kz3=diag(kz3);
% ------------------
    
[R11 R12 R21 R22]=RmatrixTM(xt,epst,d,lambda,alpha,eps1,N,M,L,nlayer);
        
% ------------------------------------------------------------------
% Finding reflection and transmission coefficients by using R matrix
% ------------------------------------------------------------------
GG3=diag((kx.^2+kz3.^2)./kz3);
GG1=diag((kx.^2+kz1.^2)./kz1);
 
 
I = eye(2*N+1);
G11=I-R22*GG1/k0;
G12=-R21*GG3/k0;
G21=-R12*GG1/k0;
G22=I-R11*GG3/k0;
G=[G11 G12;G21 G22];
 
deltai0=zeros(2*N+1,1);
deltai0(N+1)=1;
b=0;
uprim=(k1*cos(alpha)*ux-k1*sin(alpha)*uz);
b(1:2*N+1)=(-ux*I+uprim*R22/k0)*deltai0;
b(2*N+2:4*N+2)=(uprim*R12/k0)*deltai0;
RT=G\b.';
 
Rx=RT(1:2*N+1).';
Tx=RT(2*N+2:4*N+2).';
 
Rz = - kx .* Rx ./ kz1;
Tz = - kx .* Tx ./ kz3;
% ---------------------------------------------------------------
    
% --------------------------------------
% Evaluation of Diffraction efficiencies
% --------------------------------------
DE1=zeros(1,2*N+1);
DE3=zeros(1,2*N+1);
for j=-N:N
    DE1(j+N+1)=-real(kz1(j+N+1)/k1/cos(alpha))*(abs(Rx(j+N+1))^2+abs(Rz(j+N+1))^2);
    DE3(j+N+1)=real(kz3(j+N+1)/k1/cos(alpha))*(abs(Tx(j+N+1))^2+abs(Tz(j+N+1))^2);
end
dev=sum(DE3+DE1)-1;  % Energy balance
% ------------------------------------
 
PhaseR=angle(eps1*Rx(N+1)/kz1(N+1));  % N+1 ---> zeroth order
Reflection=sum(DE1);

Rmin1_0_p1=fliplr(DE1(N:N+2));

Tmin1_0_p1=fliplr(DE3(N:N+2));

%display(PhaseR)
%display(Reflection)
disp(Rmin1_0_p1)
disp(Tmin1_0_p1)



