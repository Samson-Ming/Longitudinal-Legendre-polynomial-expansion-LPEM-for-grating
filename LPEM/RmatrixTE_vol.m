function [r11 r12 r21 r22]=RmatrixTE_vol(epsx,d,lambda,alpha,eps1,N,M)

% This function calculates R matrix of a volume grating for TE polarization.
% 
% epsx: permittivity profile (Example: x=(0:100)/101; epsx = 1.5+0.02*cos (2*pi*x); )
% d: Grating thickness (Normalized to the grating period)
% lambda: Vacuum wavelength (Normalized to the grating period)
% alpha: incident angle (in radian)
% eps1: Permittivity of incident medium
% N: Truncation order
% M: Number of Legendre polynomial terms
% 
% The notation used here is very similar to following reference:
% A. Khavasi, A. K. Jahromi and K. Mehrany, "Longitudinal Legendre 
% polynomial expansion of electromagnetic fields for analysis of arbitrary 
% shaped gratings" J. Opt. Soc. Am. B, vol. 25, pp. 1564-1573 (2008).
% Equation numbers used in this code are also related to above reference

k0 = 2 * pi / lambda;
k1 = k0 * eps1 ^ 0.5;
kx = k1 * sin(alpha) - (-N:N) * 2 * pi; % see Eq. (4)

% ----------------------------------------------------
% Computing Fourier expansion of permittivity function
% ----------------------------------------------------
eps=fastft(epsx,2*N).';
eps=toeplitz(eps(2*N+1:4*N+1,1),eps(2*N+1:-1:1,1));
% ----------------------------------------------------

Q=zeros((M-1)*(2*N+1)*2,M*(2*N+1)*2);
Kx2=diag(kx.^2)/k0^2;
k0I=eye(2*N+1)*k0;
% -----------------------------------------------------------------------
% In this part the parameter "layer" is approximated for accurate results
% -----------------------------------------------------------------------
mylayer=1;
myd=d;
while myd>1e-2
    myd=myd/2;
    mylayer=mylayer+1;
end
layer=mylayer+2;

d=d/2^(layer-1);

% ------------------------------------------------------------
% Construction of "Q" matrix: see Eqs. (17)-(20) and (A1)-(A4)
% ------------------------------------------------------------
Q12=-i*k0I;
Q21=i*k0*(Kx2-eps);
for m=1:M-1
    Q(m+(0:2*N)*(M-1)+(2*N+1)*(M-1),m+(0:2*N)*M)=Q21;
    Q(m+(0:2*N)*(M-1),m+(0:2*N)*M+(2*N+1)*M)=Q12;
    for j=-N:N   
        for n=m+1+(j+N)*M:2:(j+N)*M+M
            Q(m+(j+N)*(M-1),n)=2/d*(2*m-1);
            Q(m+(j+N)*(M-1)+(2*N+1)*(M-1),n+(2*N+1)*M)=2/d*(2*m-1);
        end
    end
end
% ------------------------------------------------

% --------------------------------------------------------
% Construction of "chi" and "psi": see Eqs. (A11) and (A12)
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

% -------------------------------------
% Calculation of r matrix: see Eq. (41)
% -------------------------------------
    r=psi*([Q; chi]\[Zero; II]);
    r11=r(1:(2*N+1),1:(2*N+1));
    r12=r(1:(2*N+1),(2*N+2):(4*N+2));
    r21=r((2*N+2):(4*N+2),1:(2*N+1));
    r22=r((2*N+2):(4*N+2),(2*N+2):(4*N+2));
% -------------------------------------

% ------------------------------------------
% Successive R matrix propagation algorithm
% ------------------------------------------
for j=2:layer
    Z=inv(r22-r11);
    ZR21=Z*r21;
    ZR12=Z*r12;
    Rx11=r11-r12*ZR21;
    Rx12=r12*(ZR12);
    Rx21=-r21*ZR21;
    Rx22=r22+r21*(ZR12);
    r11=Rx11;r12=Rx12;r21=Rx21;r22=Rx22;
end
% --------------------------------------------------------