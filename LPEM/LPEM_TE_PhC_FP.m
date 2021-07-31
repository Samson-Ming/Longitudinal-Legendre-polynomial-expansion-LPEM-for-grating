% In The Name of GOD
% 
% This code has been written for analysis of conventional photonic
% crystals. The output of this program is the field profile and diffraction 
% efficiencies. If you only want to find the diffraction efficiencies, it
% is recommended to use LPEM_TE_PhC.
% The incident wave is TE polarized.
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
% 
% The outputs of this program are diffraction efficiencies which are stored
% in vectors 'DE1' and 'DE3' for reflected and transmitted waves,
% respectively. Besides, this code gives the field distribution in the 
% following matrices : 'Ey', 'Hx' and 'Hz'.

%%%%%%%
%INPUTS
%%%%%%%
close all
clear
% ----------------------
% Convergence parameters
% ----------------------
% It is recommended not to change "M = 6" and "nlayer = 5". You can achieve
% convergence by adjusting 'N' and 'L'.
N=4;     %Truncation order
L=300;        % Number of layers in a unit cell (also number of points in
             % which field profile is calculated).
Nperiod = 4; % Number periods in which field profile is plotted.
M=6;        % Number of polynomial terms in each layer.
nlayer=5;  % This parameter increase the integration accuracy

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
% arbitrary number to "radi" and "epsrod" for a homogenous layer
radi=[.2 0.2 .2]; %Rod's radius
epsrod=[1 1 1]*3.4^2;   %Rod's permittivity
epssub=[1 1 1];  %Substrate's permittivity
d_sublayer=[1 1 1];  %thickness of a sublayer of main layers (Normalized to
%                     the grating period)
dtotal = sum(d_sublayer.*N_sublayers);

% ----------------------------------------------------
% Construction of 'epst' and 'xt' for photonic crystal
% ----------------------------------------------------
j=0;
dz=0;
centerplace = 0;
r= 0;
for mj=1:length(layertype)
    if layertype(mj)==1
        dz = dz + d_sublayer(mj)*N_sublayers(mj);
    end
    if layertype(mj)==2
        for k = 1 : N_sublayers(mj)
            j=j+1;
            centerplace(1,j)=1/2;
            centerplace(2,j)= dz + d_sublayer(mj)/2;
            r(j) = radi(mj);
            dz = dz + d_sublayer(mj);
        end
    end
    if layertype(mj)==3
        for k = 1 : N_sublayers(mj)
            if mj*k == 1
                j=j+1;
                centerplace(1,j) = 0;
                centerplace(2,j) = dz;
                r(j) = radi(mj);
            end
            j=j+1;
            centerplace(1,j)=1/2;
            centerplace(2,j)= dz + d_sublayer(mj)/2;
            r(j) = radi(mj);
            dz = dz + d_sublayer(mj);
            j=j+1;
            centerplace(1,j) = 0;
            centerplace(2,j) = dz;
            r(j) = radi(mj);
        end
    end
    if layertype(mj)==4
        for k = 1 : N_sublayers(mj)
            if mj*k == 1
                j=j+1;
                centerplace(1,j) = 0;
                centerplace(2,j) = dz;
                r(j) = radi(mj-1);
            end
            j=j+1;
            centerplace(1,j)=1/2;
            centerplace(2,j)= dz + d_sublayer(mj)/2;
            r(j) = radi(mj);
            dz = dz + d_sublayer(mj);
            j=j+1;
            centerplace(1,j) = 0;
            centerplace(2,j) = dz;
            r(j) = radi(mj-1);
        end
    end
    
end
z = LPEM_zgen15(dtotal,L*nlayer);
zlen = length(z);

xt=zeros(4,zlen);
epst=zeros(5,zlen);

xt(1,:)=0;
xt(2,:)=.5;
xt(3,:)=.5;
xt(4,:)=1;

mj=1;
dendlayer =d_sublayer(mj)*N_sublayers(mj);
for k=1:zlen
    if z(k)<dendlayer
        epst(1,k)=epsrod(mj);
        epst(2,k)=epssub(mj);
        epst(3,k)=epsrod(mj);
        epst(4,k)=epssub(mj);
        epst(5,k)=epsrod(mj);
    else
        mj=mj+1;
        dendlayer =dendlayer + d_sublayer(mj)*N_sublayers(mj);
        if mj <= length(epsrod)
            epst(1,k)=epsrod(mj);
            epst(2,k)=epssub(mj);
            epst(3,k)=epsrod(mj);
            epst(4,k)=epssub(mj);
            epst(5,k)=epsrod(mj);
        end
    end
    for j=1:length(r)
        if abs(z(k)-centerplace(2,j))<r(j)
            if centerplace(1,j)==0
                xt(1,k)=sqrt(r(j)^2-(z(k)-centerplace(2,j))^2);
                xt(4,k)=1-xt(1,k);
            else
                xt(2,k)=centerplace(1,j)-sqrt(r(j)^2-(z(k)-centerplace(2,j))^2);
                xt(3,k)=1-xt(2,k);
            end
        end 
    end
end
% -----------------------------------------------------------

% --------------------------------------------------
% Main part of the code that analyzes the structure
% --------------------------------------------------
[Ey Hx Hz DE1 DE3 x zz]=FP_TE(Nperiod,xt,epst,dtotal,lambda,alpha,eps1,eps3,N,M,L,nlayer);
% ----------------------------------------------

absEy=abs(Ey);
absHx=abs(Hx);
absHz=abs(Hz);
realEy=real(Ey);
realHx=real(Hx);
realHz=real(Hz);

figure(1)
subplot(3,1,1)
imagesc (x,zz,realHz.');shading interp

sizext = size(xt); 
hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt(j,:)+k,z,'k')
    end
end
xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it H_z')

subplot(3,1,2)
imagesc (x,zz,realHx.');shading interp

hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt(j,:)+k,z,'k')
    end
end

xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it H_x')

subplot(3,1,3)
imagesc (x,zz,realEy.');shading interp

hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt(j,:)+k,z,'k')
    end
end

xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it E_y')
hold off

sum(DE3)