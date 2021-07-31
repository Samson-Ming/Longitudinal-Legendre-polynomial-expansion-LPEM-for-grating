% In The Name of GOD
% 
% This code has been written for analysis of multilayer lamellar gratings.
% The incident wave is TE polarized.
% 
% In the first section (Inputs) you should determine the structure and
% adjust convergence parameters. Here, an example of defining a structure
% is given.
% 
% Example: 
% % First layer:
% xt(:,1)=[0.25 0.75].';
% epst(:,1)=[1 4 1].';
% dl(1)=d1;
% 
% % Second layer: This layer is a homogenous layer:
% xt(:,2)=[0.25 0.75].'; % for homogenous layer it is not important how to
%                          define "xt(:,2)". The only thing that you should
%                          be careful is its size which must be compatible
%                          with other layers.
% epst(:,2)=[1 1 1].';
% dl(2)=d2;
% 
% % Third layer:
% xt(:,3)=[0.25 0.75].';
% epst(:,3)=[1 3 1].';
% dl(3)=d3;
% 
% 
% Incident plane wave(TE)
%                   \    |
%                    \   |
%                     \  |   Region 1 (eps1)
%                      \ |
%       -     ----      ----      ----      ----
%      |     |4444|    |4444|    |4444|    |4444|
%   d1<      |4444|    |4444|    |4444|    |4444|
%      |     |4444|    |4444|    |4444|    |4444|
%       -     ----      ----      ----      ----   -
%                                                   |
%                             Air                    >d2
%                                                   |
%       -     ----      ----      ----      ----   -
%      |     |3333|    |3333|    |3333|    |3333|
%   d3<      |3333|    |3333|    |3333|    |3333|
%      |     |3333|    |3333|    |3333|    |3333|
%       -     ----      ----      ----      ----
%                    <-------->            <---->
%                 Grating Period(=1)         0.5
%                  
%                   Region 3 (eps3)
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
% ----------------------
% Convergence parameters
% ----------------------
N=40;     %Truncation order
M=2;        % Number of polynomial terms.

% ---------------------------
% Global structure parameters
% ---------------------------
lambda=2.6;   %Vacuum wavelength (Normalized to the grating period,d=1)
eps1=1;     %permittivity of incident medium
eps3=1;    %permittivity of transmission medium
alpha=0*pi/180+1e-6;% Angle of incidence (in radian)

% --------------------------
% definition of grating's layers
% --------------------------
%xt(:,m) are places of transitions in m'th layer and epst(n,m) is the
%permittivity of grating between xt(n-1,m) and xt(n,m). dl is a vector whose
%m'th element is the thickness (Normalized to the grating period) of m's layer

xt(:,1)=[.25 .75].';
epst(:,1)=[1 4 1].';
%epst(:,1)=[1 4-1i 1].';
dl(1)=.5;

% This layer is a homogenous layer:
xt(:,2)=1; %???
epst(:,2)=1;
dl(2)=.2;

xt(:,3)=[.25 .75].';
epst(:,3)=[1 4 1].';
dl(3)=.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% -------------------
% initial computations
% -------------------
k0=2*pi/lambda;
k1=k0*(eps1^.5);
k3=k0*(eps3^.5);
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

Nlayer = length(dl);
for lcounter=1:Nlayer
    d=dl(lcounter);
    if max(epst(:,lcounter))==min(epst(:,lcounter))
        [r11 r12 r21 r22]=RmatrixTE_hom(epst(1,lcounter),d,lambda,alpha,eps1,N);
    else
        [r11 r12 r21 r22]=RmatrixTE_lam(xt(:,lcounter),epst(:,lcounter),d,lambda,alpha,eps1,N,M);
    end
%   -------------------------------
%   R matrix Propagation algorithm
%   -------------------------------
    if lcounter==1
         R11=r11;
         R12=r12;
         R21=r21;
         R22=r22;
    else
         Z=inv(r22-R11);
         Rx11=r11-r12*Z*r21;
         Rx12=r12*Z*R12;
         Rx21=-R21*Z*r21;
         Rx22=R22+R21*Z*R12;
         R11=Rx11;R12=Rx12;R21=Rx21;R22=Rx22;
    end
%   -------------------------------
end

I=eye(2*N+1);
G11=I+R22*Kz1/k0;
G12=R21*Kz3/k0;
G21=R12*Kz1/k0;
G22=I+R11*Kz3/k0;
G=[G11 G12;G21 G22];

deltai0=zeros(2*N+1,1);
deltai0(N+1)=1;
b=0;
uprim=(-k1*cos(alpha));
b(1:2*N+1)=(-I+uprim*R22/k0)*deltai0;
b(2*N+2:4*N+2)=(uprim*R12/k0)*deltai0;

RT=G\b.';
Ry=RT(1:2*N+1);
Ty=RT(2*N+2:4*N+2);
    
    
% --------------------------------------
% Evaluation of Diffraction efficiencies
% --------------------------------------
DE1=zeros(1,2*N+1);
DE3=zeros(1,2*N+1);
for j=-N:N
    DE1(j+N+1)=-real(kz1(j+N+1)/k1/cos(alpha))*(abs(Ry(j+N+1))^2);
    DE3(j+N+1)=real(kz3(j+N+1)/k1/cos(alpha))*(abs(Ty(j+N+1))^2);
end
dev=sum(DE3+DE1)-1; % Energy balance
% ------------------------------------
    
sum(DE1)

sum(DE3)
