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
%N=20;     %Truncation order
%L=50;       % Number of layers in a unit cell.
M=3;        % Number of polynomial terms in each layer.
nlayer = 1;  % this parameter increase the integration accuracy
% ---------------------------
% Global structure parameters
% ---------------------------
%grating period must be normalized to 1
lambda=1;     %Vacuume wavelength (Normalized to grating period)
d=1;       % Grating thickness (Normalized to grating period)
eps1=1;     %permittivity of incident medium
%eps3=1.72^2;     %permittivity of transmission medium
eps3=(1-5i)^2;
alpha=45*pi/180;% Angle of incidence (in radian)
Phi=60*pi/180;
% ------------------------
% Specifying surface type
% ------------------------
% Three different surfaces can be defined here: sinusoidal and trapezoidal
% surfaces.
%[xt,epst]=trapezoidalsurface15(0,0.25,0.75,1,eps1,eps3,d,L*nlayer);
% [xt,epst]=sinusoidalsurface15(eps1,eps3,d,L*nlayer);
%[xt,epst]=overhang_sinusoidalsurface15(eps1,eps3,d,Phi,L*nlayer);
 
RP_ref=[0.726026207305534,0.0820952701445662];
%TP_ref=[];

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =10;         % lowest number of modes
nMax_u    =100;         % highest number of modes
nMax_step = 10;          % length of the step - modes

N_l    =10;           % lower number of layers,40
N_u    =100;           % upper number of layers
N_step = 10;             % length of the step - layers

tic                   % start time count
RPvec=[];
%TSvec=[];

 for L = N_l:N_step:N_u
[xt,epst]=overhang_sinusoidalsurface15(eps1,eps3,d,Phi,L*nlayer);
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
     clear R11 R12 R21 R22 RT G b  
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

%Ref=fliplr(DE1(N+1:N+2));
%display(Ref)

% ------------------------------------
 
RP=fliplr(DE1(N+1:N+2));

%TP=[];

RPvec((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step,:)=RP;
%TPvec((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step,:)=TS;

%RT=[RP,TP];
%RT_ref=[RP_ref,TP_ref];

RT=RP;
RT_ref=RP_ref;

     %error_max((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
     error((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
                            
     c_time((N-nMax_l+nMax_step)/nMax_step,(L-N_l+N_step)/N_step) = toc;
    end
end
       tot_Run_time=sum(sum(c_time))
       
%% Save results to a file
filename = 'Metal_overhang_relief_sinusodial_h_Phi=60_TM_1percent_LPEM.mat';
%filename = 'Metal_overhang_relief_sinusodial_h_Phi=60_TM_1percent_LPEM_N=41_L=40_3_1.mat';
save(filename)       
