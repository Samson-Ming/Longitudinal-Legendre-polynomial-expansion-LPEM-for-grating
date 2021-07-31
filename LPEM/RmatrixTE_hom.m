function [R11 R12 R21 R22]=RmatrixTE_hom(epsilon,d,lambda,alpha,eps1,N)

% This function calculates R matrix of a homogenous layer for TE polarization.
% 
% epsilon: Permittivity
% d: Layer thickness (Normalized to the grating period)
% lambda: Vacuum wavelength (Normalized to the grating period)
% alpha: incident angle (in radian)
% eps1: Permittivity of incident region
% N: Truncation order

k0 = 2 * pi / lambda;
k1 = k0 * eps1 ^ 0.5;
kx = k1 * sin(alpha) - (-N:N) * 2 * pi;
k2=k0*epsilon^.5;
kz2=sqrt(k2^2-kx.^2);

mylayer=1;
myd=d;
while N*myd>10
    myd=myd/2;
    mylayer=mylayer+1;
end

Kz2=diag(kz2);
X1=diag(1./(2*i*sin(kz2*myd)));
X2X1=diag(-i./tan(kz2*myd));
X3=diag(1./(kx.^2/k0^2-epsilon));
r11=-Kz2.*X2X1.*X3/k0;
r12=2*Kz2.*X1.*X3/k0;
r21=-r12;
r22=-r11;

for j=2:mylayer
    Z=inv(r22-r11);
    ZR21=Z*r21;
    ZR12=Z*r12;
    Rx11=r11-r12*ZR21;
    Rx12=r12*(ZR12);
    Rx21=-r21*ZR21;
    Rx22=r22+r21*(ZR12);
    r11=Rx11;r12=Rx12;r21=Rx21;r22=Rx22;
end
R11=r11;
R12=r12;
R21=r21;
R22=r22;