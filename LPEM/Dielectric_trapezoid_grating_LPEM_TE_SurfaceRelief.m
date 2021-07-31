% In The Name of GOD
% 
% This code has been written for analysis of longitudinally nonhomogenous
% volume gratings. The incident wave is TE polarized.
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
% Incident plane wave(TE)
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
% ----------------------
% Convergence parameters
% ----------------------
% It is recommended not to change "M" and "nlayer".
%N=30;      %Truncation order
%L=5;       % Number of layers in a unit cell.
M=7;        % Number of polynomial terms in each layer.
nlayer = 5;  % this parameter increase the integration accuracy
% ---------------------------
% Global structure parameters
% ---------------------------
%grating period must be normalized to 1
lambda=1/1.5;     %Vacuum wavelength (Normalized to the grating period)
d=1;      % Grating thickness (Normalized to the grating period)
eps1=1;     %Permittivity of incident medium
eps3=1.5^2;     %Permittivity of transmission medium
%eps3=(1-5i)^2;     %Permittivity of transmission medium
%alpha=15*pi/180;% Angle of incidence (in radian)
alpha=0*pi/180;% Angle of incidence (in radian)
% ------------------------
% Specifying surface type
% ------------------------
% Two different surfaces can be defined here: sinusoidal and trapezoidal
% surfaces.
%Phi=0*(pi/180);
%Phi=15*(pi/180);
%[xt,epst]=trapezoidalsurface15(1/4,5/8,7/8,3/4,eps1,eps3,d,L*nlayer);
% [xt,epst]=sinusoidalsurface15(eps1,eps3,d,L*nlayer);
 
RS_ref=[0.0116110339181890	0.00397254128892009	0.00269456222006578];
TS_ref=[0.0984302860119821	0.476629669959978	0.0994663070054246	0.268814629169694	0.0384145463940903];

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =7;         % lowest number of modes
nMax_u    =7;         % highest number of modes
nMax_step = 1;          % length of the step - modes

N_l    =2;           % lower number of layers
N_u    = 2;           % upper number of layers
N_step = 1;             % length of the step - layers

tic                   % start time count
RSvec=[];
TSvec=[];

RSvec=[];
TSvec=[];
 for L = N_l:N_step:N_u
[xt,epst]=trapezoidalsurface15(1/4,5/8,7/8,3/4,eps1,eps3,d,L*nlayer);
for N =  nMax_l:nMax_step:nMax_u
    tic
    N
    L
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
    
[R11 R12 R21 R22]=RmatrixTE(xt,epst,d,lambda,alpha,eps1,N,M,L,nlayer);
        
% ---------------------------------------------------------------
% Finding reflection and transmission coefficients using R matrix
% ---------------------------------------------------------------
I = eye(2*N+1);
G11=I+R22*Kz1/k0;
G12=R21*Kz3/k0;
G21=R12*Kz1/k0;
G22=I+R11*Kz3/k0;
G=[G11 G12;G21 G22];

deltai0=zeros(2*N+1,1);
deltai0(N+1)=1;
uprim=(-k1*cos(alpha));
b(1:2*N+1)=(-I+uprim*R22/k0)*deltai0;
b(2*N+2:4*N+2)=(uprim*R12/k0)*deltai0;


RT=G\b.';
Ry=RT(1:2*N+1);
Ty=RT(2*N+2:4*N+2);
% ---------------------------------------------------------------
    clear R11 R12 R21 R22 RT G b
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
 
%PhaseR=angle(eps1*Rx(N+1)/kz1(N+1));  % N+1 ---> zeroth order
%Reflection=sum(DE1);

RS=fliplr(DE1(N:N+2));

TS=fliplr(DE3(N-1:N+3));

RSvec((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step,:)=RS;
TSvec((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step,:)=TS;

RT=[RS,TS];
RT_ref=[RS_ref,TS_ref];

     %error_max((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
     error((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
                            
     c_time((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step) = toc;
    end
end
       tot_Run_time=sum(sum(c_time))
       
%% Save results to a file
%filename = 'dielectric_overhang_relief_trapezoid_TE_1percent_LPEM.mat';
filename = 'dielectric_overhang_relief_trapezoid_TE_1percent_LPEM_N=15_L=2.mat';
save(filename)       

