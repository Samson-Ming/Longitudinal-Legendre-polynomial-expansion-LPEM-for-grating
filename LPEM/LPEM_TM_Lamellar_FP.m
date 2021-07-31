% In The Name of GOD
% 
% This code has been written for analysis of multilayer lamellar gratings.
% The output of this program is the field profile and diffraction 
% efficiencies. If you only want to find the diffraction efficiencies, it
% is recommended to use LPEM_TM_Lamellar. The incident wave is TM polarized.
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
% Incident plane wave(TM)
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
% respectively. Besides, this code gives the field distribution in the 
% following matrices: 'Hy', 'Ex' and 'Ez'.
 
%%%%%%%
%INPUTS
%%%%%%%
clear
close all
% ----------------------
% Convergence parameters
% ----------------------
N=10;     %Truncation order
L=100;        % Number of layers in a unit cell (also number of points in
             % which field profile is calculated).
Nperiod = 3; % Number periods in which field profile is plotted.
M=3;        % Number of polynomial terms in each layer.
 
% ---------------------------
% Global structure parameters
% ---------------------------
lambda=2.6;   %Vacuum wavelength (Normalized to the grating period)
eps1=1;     %Permittivity of incident medium
eps3=1;    %Permittivity of transmission medium
alpha=0*pi/180+1e-6;% Angle of incidence (in radian)
 
% ------------------------------
% definition of grating's layers
% ------------------------------
% xt(:,m) are places of transitions in m'th layer and epst(n,m) is the
% permittivity of grating between xt(n-1,m) and xt(n,m). dl is a vector whose
% m'th element is the thickness (Normalized to grating period) of m's
% layer.
 
xt(:,1)=[.25 .75].';
epst(:,1)=[1 4 1].';
dl(1)=0.5;
 
% This layer is a homogenous layer:
xt(:,2)=0.5;
epst(:,2)=1;
dl(2)=0.2;
% 
xt(:,3)=[.25 0.75].';
epst(:,3)=[1 4 1].';
dl(3)=0.5;
%------------------------------------
 
nlayer = 1;
 
z = LPEM_zgen15(sum(dl),L*nlayer);
zlen = length(z);
 
sizeepst=size(epst);
sizext=size(xt);
 
xt2=zeros(sizext(1),zlen);
epst2=zeros(sizeepst(1),zlen);
 
 
mj=1;
dendlayer =dl(mj);
for k=1:zlen
    if z(k)<dendlayer
        epst2(:,k)=epst(:,mj);
        xt2(:,k)=xt(:,mj);
    else
        mj=mj+1;
        dendlayer =dendlayer + dl(mj);
        if mj <= length(dl)
            epst2(:,k)=epst(:,mj);
            xt2(:,k)=xt(:,mj);
        end
    end
end
 
 
 
d=sum(dl);
 
% tic;
% [Hy Ex Ez DE1 DE3 x zz]=FP_TM(Nperiod,xt2,epst2,d,lambda,alpha,eps1,eps3,N,M,L,nlayer);
% sum(DE1)
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
Kx = diag(kx);
% ------------------
 
R11FP=zeros(2*N+1,2*N+1,L);
R12FP=zeros(2*N+1,2*N+1,L);
R21FP=zeros(2*N+1,2*N+1,L);
R22FP=zeros(2*N+1,2*N+1,L);
r21FP=zeros(2*N+1,2*N+1,L);
ZFP=zeros(2*N+1,2*N+1,L);
 
dfir = 0;
for mj = 1 : length(dl)
    LL(mj) = floor(dl(mj)*L/d);
    dz = dl(mj) / LL(mj);
    
    if mj == 1
        zz = (0 : dz : dl(1) - dz);
    else
        zz = [zz dfir:dz:dfir+dl(mj)-dz];
    end
    dfir = dfir + dl(mj);
end
L = sum(LL);
zz (L+1) = d;
% return
Lend = 1;
for mj = 1 : length(dl)
    [r11 r12 r21 r22]=RmatrixTM_lam(xt(:,mj),epst(:,mj),dl(mj)/LL(mj),lambda, ... 
                                alpha,eps1,N,M);
    
for lcounter = Lend : Lend+LL(mj)-1
    if (lcounter==1)
        R11=r11;
        R12=r12;
        R21=r21;
        R22=r22;
        Z=eye(size(r22));
    else
        Z=inv(r22-R11);
        Rx11=r11-r12*Z*r21;
        Rx12=r12*Z*R12;
        Rx21=-R21*Z*r21;
        Rx22=R22+R21*Z*R12;
        R11=Rx11;R12=Rx12;R21=Rx21;R22=Rx22;
    end
    R11FP(:,:,lcounter)=R11;
    R12FP(:,:,lcounter)=R12;
    R21FP(:,:,lcounter)=R21;
    R22FP(:,:,lcounter)=R22;
    ZFP(:,:,lcounter)=Z;
    r21FP(:,:,lcounter)=r21;
end
Lend = Lend + LL(mj);
end
% ---------------------------------------------------------------
% Finding reflection and transmission coefficients using R matrix
% ---------------------------------------------------------------
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
% ------------------------------------
 
Uy=zeros(2*N+1,L);
Sx=Uy;
Sz=Uy;
 
Sx(:,1)=Rx;
Sx(N+1,1)=Sx(N+1,1)+ux;
 
Uy(:,1)=kz1.*Rx-kx.*Rz;
Uy(N+1,1)=Uy(N+1,1)+k1*cos(alpha)*ux-kx(N+1)*uz;
Uy(:,1)=Uy(:,1)/k0;
 
Uyd=(kz3.*Tx-kx.*Tz)/k0;
Sxd=Tx;
Uy(:,L+1)=Uyd;
Sx(:,L+1)=Sxd;
 
 
for j=L+1:-1:3
    R11=R11FP(:,:,j-2);
    R12=R12FP(:,:,j-2);
    Z=ZFP(:,:,j-1);
    r21=r21FP(:,:,j-1);
    Uy(:,j-1)=Z*(-r21*Uy(:,j)+R12*Uy(:,1));
    Sx(:,j-1)=R11*Uy(:,j-1)+R12*Uy(:,1);
end
 
mj=1;
dendlayer =dl(mj);
for j = 1:L+1
    if zz(j)<=dendlayer
        epsmn=(fastftptrain(xt(:,mj),epst(:,mj),2*N)).';
    else
        mj = mj + 1;
        dendlayer = dendlayer + dl(mj);
        epsmn=(fastftptrain(xt(:,mj),epst(:,mj),2*N)).';
        
    end
    eps_M=toeplitz(epsmn(2*N+1:4*N+1),epsmn(2*N+1:-1:1));
    Sz(:,j) = -inv(eps_M)*Kx*Uy(:,j)/k0;
end
 
xlen=(L+1)*Nperiod;
x = (0:xlen)*Nperiod/(xlen+1);
Hy = zeros(xlen+1,L+1);
Ex = zeros(xlen+1,L+1);
Ez = zeros(xlen+1,L+1);
eta0 = sqrt(4*pi*1e-7/(8.85e-12));
 
 
for j=-N:N
      Ex=Ex+(exp(-i*kx(j+N+1)*x).')*Sx(j+N+1,:);
      Hy=Hy+(exp(-i*kx(j+N+1)*x).')*Uy(j+N+1,:)/eta0;
      Ez=Ez+(exp(-i*kx(j+N+1)*x).')*Sz(j+N+1,:);
end
% toc;
sum(DE1)
 
absHy=abs(Hy);
absEx=abs(Ex);
absEz=abs(Ez);
realHy=real(Hy);
realEx=real(Ex);
realEz=real(Ez);
 
figure(1)
subplot(3,1,1)
imagesc (x,zz,realEz.');shading interp
% colormap(gray)
sizext = size(xt); 
hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt2(j,:)+k,z,'k')
    end
end
xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it E_z')
 
subplot(3,1,2)
imagesc (x,zz,realEx.');shading interp
% colormap(gray)
hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt2(j,:)+k,z,'k')
    end
end
 
xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it E_x')
 
subplot(3,1,3)
imagesc (x,zz,realHy.');shading interp
% colormap(gray)
% quiver(absEz,absEx)
hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt2(j,:)+k,z,'k')
    end
end
 
xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it H_y')
hold off
% figure(2)
% quiver(x,z,powerx.',powerz.');shading interp