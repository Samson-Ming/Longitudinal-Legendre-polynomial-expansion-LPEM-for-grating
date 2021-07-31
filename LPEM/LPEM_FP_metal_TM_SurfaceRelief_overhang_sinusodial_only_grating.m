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
% ----------------------
% Convergence parameters
% ----------------------
% It is recommended not to change "M" and "nlayer".
N=20;     %Truncation order
L=40;       % Number of layers in a unit cell.
M=3;        % Number of polynomial terms in each layer.
nlayer = 1;  % this parameter increase the integration accuracy
Nperiod=1;
% ---------------------------
% Global structure parameters
% ---------------------------
%grating period must be normalized to 1
abs_period=1;
lambda=1;     %Vacuume wavelength (Normalized to grating period)
d2=0.5;       % Grating thickness (Normalized to grating period)
eps1=1;     %permittivity of incident medium
%eps3=1.72^2;     %permittivity of transmission medium
eps3=(1-5i)^2;
alpha=45*pi/180;% Angle of incidence (in radian)
Phi=60*pi/180;
% ------------------------
% Specifying surface type
% ------------------------
% Two different surfaces can be defined here: sinusoidal and trapezoidal
% surfaces.
%[xt,epst]=trapezoidalsurface15(0,0.25,0.75,1,eps1,eps3,d,L*nlayer);
% [xt,epst]=sinusoidalsurface15(eps1,eps3,d,L*nlayer);
%[xt,epst,xt_original,z]=overhang_sinusoidalsurface15(eps1,eps3,d,Phi,L*nlayer);
%z=LPEM_zgen15(d,L*nlayer);

xt1=[.25 .75].';
epst1=eps1*[1 1 1].';
%dl(1)=0.1*d2;
dl(1)=0*d2;
 
% This layer is a homogenous layer:
%xt(:,2)=0.5;
%epst(:,2)=1;
dl(2)=d2;
% 
xt3=[.25 0.75].';
epst3=eps3*[1 1 1].';
%dl(3)=0.1*d2;
dl(3)=0*d2;
%------------------------------------
z = LPEM_zgen15(d2,L*nlayer);
zlen = length(z);

xt2=zeros(2,zlen);
epst2=zeros(3,zlen);

%mj=1;
%dendlayer =dl(mj);
xt_original=zeros(2,1);
z2=0;
k2=0;
for k=1:zlen
        k2=k2+1;
        %[xt,epst,xt_original,z]=overhang_sinusoidalsurface15(eps1,eps3,d,Phi,L*nlayer);
        epst2(1,k)= eps1;
        epst2(2,k) = eps3;
        epst2(3,k) = eps1;
        xt2(1,k)=acos(-(2*((dl(1)+dl(2))-z(k))/d2-1))/(2*pi/cos(Phi));
        xt2(2,k)=cos(Phi)-xt2(1,k);
    
        xt2(1,k)=sec(Phi)*xt2(1,k)+tan(Phi)*((dl(1)+dl(2))-z(k));
        xt2(2,k)=sec(Phi)*xt2(2,k)+tan(Phi)*((dl(1)+dl(2))-z(k));
    
       xt_original(1,k2)=xt2(1,k);
       xt_original(2,k2)=xt2(2,k);
       z2(k2)=z(k);
    
    %%{
    %%{
    x1=mod(xt2(1,k),1);
    x2=mod(xt2(2,k),1);
    if x1>x2
        xt2(1,k)=x2;
        xt2(2,k)=x1;
        epst2(1,k)= eps3;
        epst2(2,k) = eps1;
        epst2(3,k) = eps3;
    else
        xt2(1,k)=x1;
        xt2(2,k)=x2;
    end
%}
end

%------------------------------------
%{
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
Tr=fliplr(DE3(N:N+2));
PhaseR=angle(eps1*Rx(N+1)/kz1(N+1));  % N+1 ---> zeroth order
Reflection=sum(DE1);
 
%display(PhaseR)
%display(Reflection)
display(Tr)
%}
% --------------------------------------------------
% Main part of the code that analyzes the structure
% --------------------------------------------------
[Hy Ex Ez DE1 DE3 x zz]=FP_TM(Nperiod,xt2,epst2,d2,lambda,alpha,eps1,eps3,N,M,L,nlayer);
% --------------------------------------------------
absHy=abs(Hy);
absH=absHy;
absEx=abs(Ex);
absEz=abs(Ez);
%absE=sqrt(absEx.^2+absEz.^2);
realHy=real(Hy);
realEx=real(Ex);
realEz=real(Ez);

%{
figure(1)
subplot(3,1,1)
imagesc (x,zz,realEz.');shading interp

sizext = size(xt); 
hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt(j,:)+k,z,'k')
    end
end
xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it E_z')

subplot(3,1,2)
imagesc (x,zz,realEx.');shading interp

hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot(xt(j,:)+k,z,'k')
    end
end

xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it E_x')

subplot(3,1,3)
imagesc (x,zz,realHy.');shading interp

hold on 
for j = 1:sizext(1)
    for k=0:Nperiod-1
        plot((xt(j,:)+k)*abs_period,(d-z)*abs_period,'k')
    end
end

xlabel('x/\Lambda_G')
ylabel('z/\Lambda_G')
title('\it H_y')
hold off
%}

figure,axis equal
%imagesc (x*abs_period,((dl(1)+dl(2))-zz)*abs_period,absH.');
pcolor (x*abs_period,((dl(1)+dl(2))-zz)*abs_period,120*pi*absH.');
shading interp
colormap jet
colorbar

sizext = size(xt2); 
hold on 
%%{
for k=0:Nperiod-1
    for j = 1:sizext(1)
        plot((xt_original(j,:)+k-1)*abs_period,((dl(1)+dl(2))-z2)*abs_period,'w-','LineWidth',1.5); 
        plot((xt_original(j,:)+k)*abs_period,((dl(1)+dl(2))-z2)*abs_period,'w-','LineWidth',1.5);
        plot((xt_original(j,:)+k+1)*abs_period,((dl(1)+dl(2))-z2)*abs_period,'w-','LineWidth',1.5); 
    end
    plot([(xt_original(2,1)+k-1),(xt_original(1,1)+k-1)]*abs_period,[((dl(1)+dl(2))-z2(1)),((dl(1)+dl(2))-z2(1))]*abs_period,'w-','LineWidth',1.5); 
    plot([(xt_original(2,end)+k-1),(xt_original(1,end)+k)]*abs_period,[((dl(1)+dl(2))-z2(end)),((dl(1)+dl(2))-z2(end))]*abs_period,'w-','LineWidth',1.5);
    
end
%}
%{
xxt=[xt_original(1,:),fliplr(xt_original(2,:))];
yyt=[((dl(1)+dl(2))-z2),fliplr(((dl(1)+dl(2))-z2))];
for k=0:Nperiod-1
    plot((xxt+k-1)*abs_period,yyt*abs_period,'k-','LineWidth',1.5); 
    plot((xxt+k)*abs_period,yyt*abs_period,'k-','LineWidth',1.5); 
    plot((xxt+k+1)*abs_period,yyt*abs_period,'k-','LineWidth',1.5); 
end
%}
xlabel('x')
ylabel('z')
title('\it |H|')
axis xy
hold off

dev=sum(DE3+DE1)-1; % Energy balance
% ------------------------------------
Tr=fliplr(DE3(N:N+2));
Rr=fliplr(DE1(N:N+2));

display(Tr)
display(Rr)
filename = 'Near_field_only_grating_Metal_overhang_relief_sinusodial_0.5h_Phi=60_TM_1percent_LPEM_N=41_L=40_3_1.mat';
save(filename)     