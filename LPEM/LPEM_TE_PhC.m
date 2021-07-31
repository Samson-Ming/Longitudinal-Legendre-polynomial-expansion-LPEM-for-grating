% In The Name of GOD
% 
% This code has been written for analysis of conventional photonic
% crystals. The incident wave is TE polarized.
% 
% In the first section (Inputs) you should determine the structure and
% adjust convergence parameters. Here, an example of defining a structure
% is given.
% 
% Example: A square lattice with circular rods and a defect line
% 
% layertype=[2 1 2]; %1--> homogenous, 2--> square lattice
% N_sublayers=[2 1 2]; % Number of sublayers that have same structures
% radi=[r 0 r]; %Radius of rods
% epsrod=[eps(GaAs) 1 eps(Silicon)];   %Rod's permittivity
% epssub=[1 1 1];  %Substrate's permittivity
% d_sublayer=[d1 d2 d3]; %thickness of a sublayer of main layers 
% 
% 
% Incident plane wave(TE)
%                   \    |
%                    \   |
%                     \  |   Region 1 (eps1)
%                      \ |
%        -     ------------------------------   - 
%       |          --        --        --        |
%       |        /****\    /GaAs\    /****\      |
%   d1 <         \****/    \****/    \****/      |
%       |          --        --        --        |  First main layer: a 
%        -                  Air                   > square lattice with two
%                  --        --        --        |  sublayers
%                /****\    /****\    /****\      |  
%                \****/    \****/    \****/      |
%                  --        --        --        |
%        -     ------------------------------   -
%       |                                        | 
%   d2 <                    Air                   >Second main layer: a 
%       |                                        | homogenous layer
%        -     ------------------------------   -
%       |          --        --        --        |
%       |        /****\    /*Si*\    /****\      |
%   d3 <         \****/    \****/    \****/      |
%       |          --        --        --        | First main layer: a 
%        -                  Air                   >square lattice with two
%                  --        --        --        | sublayers
%                /****\    /****\    /****\      |
%                \****/    \****/    \****/      |
%                  --        --        --        |
%              ------------------------------   -
%                        <-------->
%                     Grating Period(=1)
%                  
%                       Region 3 (eps3)
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
% It is recommended not to change "M = 6" and "nlayer = 5". You can achieve
% convergence by adjusting 'N' and 'L'.
N=4;     %Truncation order
L=20;        % Number of layers in a unit cell.
M=6;        % Number of polynomial terms in each layer.
nlayer=5;  % this parameter increase the integration accuracy
 
% ---------------------------
% Global structure parameters
% ---------------------------
lambda=3.271;     %Vacuum wavelength (Normalized to the grating period)
eps1=3.4^2;     %permittivity of incident medium
eps3=3.4^2;     %permittivity of transmission medium
alpha=2*pi/180+1e-9;% Angle of incidence(in radian)
% ------------------------------------------
% defining different main layers' parameters
% ------------------------------------------
% Notice that each layer can be consisted of several SAME sublayers. For
% example a square lattice layer can be a slice of a 2D photonic crystal
% with 3,4,... layers with same unit cells.
layertype=[2 1 2];  %1--> homogenous, 2--> square lattice, 3--> triangular lattice
                  %4--> triangular lattice defect
                   
N_sublayers=[3 1 3];
% For a homogeneous layer one must only define epssub but the length of 
% other vectors must be the same as epsub. So you can attribute an
% aribtrary number to "radi" and "epsrod" for a homogenous layer
radi=[.2 0 .2]; % Rod's radius
epsrod=[4 4 4];   % Rod's permittivity
epssub=[1 1 1]*3.4^2;  % Substrate's permittivity
d_sublayer=[1 1 1];  % thickness of a sublayer of main layers (Normalized to 
                     % the grating period)
dtotal = sum(d_sublayer.*N_sublayers);
 
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
 
 
for mj=1:length(layertype)
% ---------------------------------------------
% Evaluation of r matrix of different sublayers
% ---------------------------------------------
    n=L*nlayer;
    d=d_sublayer(mj);
    if layertype(mj)==1
        [R11 R12 R21 R22]=RmatrixTE_hom(epssub(mj),d,lambda,alpha,eps1,N);
    else
              
        if layertype(mj)==2
            [xt,epst]=SquareLatticeCircularRod15(d,radi(mj),epsrod(mj),epssub(mj),n);
        end
        if layertype(mj)==3
            [xt,epst]=TriangularLatticeCircularRod15(d,radi(mj),epsrod(mj),epssub(mj),n);
        end
        if layertype(mj)==4
            [xt,epst]=TriangularLatticeCircularRodDefect15(d,radi(mj-1),radi(mj),epsrod(mj-1),epssub(mj-1),epsrod(mj),epssub(mj),n);
        end
        
        [R11 R12 R21 R22]=RmatrixTE(xt,epst,d,lambda,alpha,eps1,N,M,L,nlayer);
        
        
    end
%   -------------------------------------------
 
    
%   ----------------------------------------------------------------
%   R matrix propagation algorithm for a main layer that consists of
%   several same sublayers.
%   ----------------------------------------------------------------
    
    r11=R11;
    r12=R12;
    r21=R21;
    r22=R22;
    for j=2:N_sublayers(mj)
        Z=inv(r22-R11);
        Rx11=r11-r12*Z*r21;
        Rx12=r12*Z*R12;
        Rx21=-R21*Z*r21;
        Rx22=R22+R21*Z*R12;
        R11=Rx11;R12=Rx12;R21=Rx21;R22=Rx22;
    end
    r11=R11;
    r12=R12;
    r21=R21;
    r22=R22;
%   ----------------------------------------------------------------
 
%   -------------------------------
%   R matrix Propagation algorithm
%   -------------------------------
    if mj==1
         RT11=r11;
         RT12=r12;
         RT21=r21;
         RT22=r22;
    else
         Z=inv(r22-RT11);
         Rx11=r11-r12*Z*r21;
         Rx12=r12*Z*RT12;
         Rx21=-RT21*Z*r21;
         Rx22=RT22+RT21*Z*RT12;
         RT11=Rx11;RT12=Rx12;RT21=Rx21;RT22=Rx22;
    end
%   -------------------------------
end
 
% -------------------------------------------
% R matrix of the structure is now calculated
% -------------------------------------------
R11=RT11;
R12=RT12;
R21=RT21;
R22=RT22;
% -------------------------------------------
% [DE1 DE3]=DETE(R11,R12,R21,R22,lambda,alpha,eps1,eps3);
    
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
 
PhaseR = angle(Ry(N+1));
Reflection = sum(DE1);
display(PhaseR)
display(Reflection)

