function [Ey Hx Hz DE1 DE3 x z]=FP_TE(Nperiod,xt,epst,d,lambda,alpha,eps1,eps3,N,M,L,nlayer)
 
% This function calculates field profile (Ey Hx Hz) and diffraction 
% efficiencies (DE1 DE3) of a grating with discontinuous profile for TE
% polarization.
% 
% 
% Nperiod : Number of periods in which field profile is calculated
% xt: Place of transition of permittivity profile
% epst: Permittivity between transitions
% 
% To understand how you can define a grating by determining 'xt' and
% 'epst', look at these functions: sinusoidalsurface15.m,
% SquareLatticeCircularRod15.m, trapezoidalsurface15.m, ...
% Notice that: 1. Number of columns of 'epst' and 'xt' must be L*nlayer
%              2. Number of rows of 'epst' must be number of rows of 'xt' + 1
%              3. The direction of coordinates is like this:
% 
%                       Region 1
%                  ------------------ z=0   |-----> x
%                    Grating Region         |     
%                  ------------------ z=d   V
%                       Region 3            z
% 
% d: Grating thickness (Normalized to grating period)
% lambda: Vacuum wavelength (Normalized to grating period)
% alpha: incident angle (in radian)
% eps1: Permittivity of incident region
% eps3: Permittivity of transmission region
% N: Truncation order
% M: Number of Legendre polynomial terms
% L: Number of layers (Also number of points in z and x direction --> 
%    size(Ey)=(L+1)*(L+1))
% nlayer: Number of Gauss-Legendre integration algorithm used in each layer
% 
% The notation used here is very similar to following reference:
% A. Khavasi, A. K. Jahromi and K. Mehrany, "Longitudinal Legendre 
% polynomial expansion of electromagnetic fields for analysis of arbitrary 
% shaped gratings" J. Opt. Soc. Am. B, vol. 25, pp. 1564-1573 (2008).
% Equation numbers used in this code are also related to above reference
 
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
 
eps=zeros(4*N+1,15*nlayer);
w=zeros(1,15*nlayer);
 
xeps=zeros(4*N+1,M,M);
xepsM=zeros(2*N+1,2*N+1,M,M);
 
% Weights of Gauss-Legendre integration algorithm
w1=0.030753241996117;
w2=0.070366047488108;
w3=0.107159220467172;
w4=0.139570677926154;
w5=0.166269205816994;
w6=0.186161000115562;
w7=0.198431485327111;
w8=0.202578241925561;
 
Kx=diag(kx);
delta=1/nlayer;
d=d/L;
P=leg_gen15(M-1,nlayer); %Legendre polynomials
 
R11FP=zeros(2*N+1,2*N+1,L);
R12FP=zeros(2*N+1,2*N+1,L);
R21FP=zeros(2*N+1,2*N+1,L);
R22FP=zeros(2*N+1,2*N+1,L);
r21FP=zeros(2*N+1,2*N+1,L);
ZFP=zeros(2*N+1,2*N+1,L);
 
I=eye(2*N+1);
K0I=-i*k0*I*2;
Kx2=i*Kx^2*2/k0;
Q=zeros((M-1)*(2*N+1)*2,M*(2*N+1)*2);
% --------------------------------------------------------
% Construction of "chi" and "psi": see Eq. (A11) and (A12)
% --------------------------------------------------------
chi=zeros(2*(2*N+1),2*(2*N+1)*M); 
psi=zeros(2*(2*N+1),2*(2*N+1)*M);
chi22=zeros((2*N+1),(2*N+1)*M);
chi12=chi22;
for j=-N:N
    chi22(j+N+1,(1+(j+N)*M):((j+N+1)*M))=(-1).^(0:(M-1));
    chi12(j+N+1,(1+(j+N)*M):((j+N+1)*M))=1;
end
chi(1:(2*N+1),((2*N+1)*M+1):2*(2*N+1)*M)=chi12;
chi((2*N+2):(4*N+2),((2*N+1)*M+1):2*(2*N+1)*M)=chi22;
psi(1:(2*N+1),1:(2*N+1)*M)=chi12;
psi((2*N+2):(4*N+2),1:(2*N+1)*M)=chi22;
% --------------------------------------------------------
Zero=zeros(2*(2*N+1)*(M-1),2*(2*N+1));
II=eye(2*(2*N+1),2*(2*N+1));
 
wx=[w1 w2 w3 w4 w5 w6 w7 w8 w7 w6 w5 w4 w3 w2 w1]*delta;
for k=1:nlayer
    w((1+15*(k-1)):15*k)=wx;
end
for lcounter=1:L
    k=(lcounter-1)*15*nlayer+1;
    
    for kk=1:15*nlayer            
%         ----------------------------------------------------
%         Computing Fourier expansion of permittivity function
%         ----------------------------------------------------
        eps(:,kk)=-i*k0*(fastftptrain(xt(:,k+kk-1),epst(:,k+kk-1),2*N)).';
%         ----------------------------------------------------
    end
% --------------------------
% Gauss-Legendre integration
% --------------------------
    for m1=1:M
        for m2=1:m1
            Pm1m2=P(m1,:).*P(m2,:).*w;
            xeps(:,m1,m2)=eps*Pm1m2.';
            xepsM(:,:,m1,m2)=toeplitz(xeps(2*N+1:4*N+1,m1,m2),xeps(2*N+1:-1:1,m1,m2));
            xepsM(:,:,m2,m1)=xepsM(:,:,m1,m2);
        end
    end
% ---------------------------
 
% -----------------------------------------------
% Construction of "Q" matrix: see Eqs. (A1)-(A4)
% -----------------------------------------------
    for m1=1:M-1
        for m2=1:M
%           Eq. (A2)
            Q(m1+(0:(M-1):2*N*(M-1))+(2*N+1)*(M-1),m2+(0:M:2*N*M))=xepsM(:,:,m1,m2);
            
        end
    end
    for m=1:M-1
%       Eq. (A1)
        Q(m+(0:2*N)*(M-1),m+(0:2*N)*M+(2*N+1)*M)=K0I/(2*m-1);
%       Next part of Eq. (A2)
        Q(m+(0:2*N)*(M-1)+(2*N+1)*(M-1),m+(0:2*N)*M)=Q(m+(0:2*N)*(M-1)+(2*N+1)*(M-1),m+(0:2*N)*M)+Kx2/(2*m-1);
        for j=-N:N
            for n=m+1+(j+N)*M:2:(j+N)*M+M
%               Eq. (A3)
                Q(m+(j+N)*(M-1),n)=4/d;
%               Eq. (A4)
                Q(m+(j+N)*(M-1)+(2*N+1)*(M-1),n+(2*N+1)*M)=4/d;
            end
        end
    end
 
% -------------------------------------
% Calculation of r matrix: see Eq. (41)
% -------------------------------------
    r=psi*([Q; chi]\[Zero; II]);
    r11=r(1:(2*N+1),1:(2*N+1));
    r12=r(1:(2*N+1),(2*N+2):(4*N+2));
    r21=r((2*N+2):(4*N+2),1:(2*N+1));
    r22=r((2*N+2):(4*N+2),(2*N+2):(4*N+2));
% -------------------------------------
 
% ------------------------------------------------------
% R matrix propagation algorithm: see Eqs. (42) and (43)
% ------------------------------------------------------
    if lcounter==1
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
% -------------------------------------------------------
        R11FP(:,:,lcounter)=R11;
        R12FP(:,:,lcounter)=R12;
        R21FP(:,:,lcounter)=R21;
        R22FP(:,:,lcounter)=R22;
        ZFP(:,:,lcounter)=Z;
        r21FP(:,:,lcounter)=r21;
end
 
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
% ------------------------------------
 
 
Ux=zeros(2*N+1,L);
Sy=Ux;
Uz=Ux;
 
Sy(:,1)=Ry;
Sy(N+1,1)=Sy(N+1,1)+1;
Uz(:,1)=(kx.').*Sy(:,1)/k0;
 
Ux(:,1)=-(kz1.').*Ry;
Ux(N+1,1)=Ux(N+1,1)-k1*cos(alpha);
Ux(:,1)=Ux(:,1)/k0;
 
Uxd=-(kz3.').*Ty/k0;
Syd=Ty;
Ux(:,L+1)=Uxd;
Sy(:,L+1)=Syd;
Uz(:,L+1)=(kx.').*Sy(:,L+1)/k0;
 
for j=L+1:-1:3
    R11=R11FP(:,:,j-2);
    R12=R12FP(:,:,j-2);
    Z=ZFP(:,:,j-1);
    r21=r21FP(:,:,j-1);
    Ux(:,j-1)=Z*(-r21*Ux(:,j)+R12*Ux(:,1));
    Sy(:,j-1)=R11*Ux(:,j-1)+R12*Ux(:,1);
    Uz(:,j-1)=(kx.').*Sy(:,j-1)/k0;
end
 
xlen=(L+1)*Nperiod;
x = (0:xlen)*Nperiod/(xlen+1);
z = linspace(0,d*L,L+1);
Ey = zeros(xlen+1,L+1);
Hx = zeros(xlen+1,L+1);
Hz = zeros(xlen+1,L+1);
eta0 = sqrt(4*pi*1e-7/(8.85e-12));
for j=-N:N
    %{
       for xx=1:xlen+1
            Ey(xx,:)=Ey(xx,:)+Sy(j+N+1,:)*exp(-i*kx(j+N+1)*x(xx));
            Hx(xx,:)=Hx(xx,:)+Ux(j+N+1,:)*exp(-i*kx(j+N+1)*x(xx))/eta0;
            Hz(xx,:)=Hz(xx,:)+Uz(j+N+1,:)*exp(-i*kx(j+N+1)*x(xx))/eta0;
       end
    %}
    Ey=Ey+(exp(-i*kx(j+N+1)*x).')*Sy(j+N+1,:);
    Hx=Hx+(exp(-i*kx(j+N+1)*x).')*Ux(j+N+1,:)/eta0;
    Hz=Hz+(exp(-i*kx(j+N+1)*x).')*Uz(j+N+1,:)/eta0;
end
