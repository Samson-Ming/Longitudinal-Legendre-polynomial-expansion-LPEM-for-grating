function [Hy Ex Ez DE1 DE3 x z]=FP_TM(Nperiod,xt,epst,d,lambda,alpha,eps1,eps3,N,M,L,nlayer)
 
% This function calculates field profile (Hy Ex Ez) and diffraction 
% efficiencies (DE1 DE3) of a grating with discontinuous profile for TM
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
% d: Grating thickness (Normalized to the grating period)
% lambda: Vacuum wavelength (Normalized to the grating period)
% alpha: incident angle (in radian)
% eps1: Permittivity of incident region
% eps3: Permittivity of transmission region
% N: Truncation order
% M: Number of Legendre polynomial terms
% L: Number of layers (Also number of points in z and x direction --> 
%    size(Hy)=(L+1)*(L+1))
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
 
 
z=LPEM_zgen15(d,L*nlayer);
 
Q11=zeros(2*N+1,2*N+1,M,M);
Q12=zeros(2*N+1,2*N+1,M,M);
Q21=zeros(2*N+1,2*N+1,M,M);
Q22=zeros(2*N+1,2*N+1,M,M);
 
 
% Weights of Gauss-Legendre integration algorithm
w1=0.030753241996117;
w2=0.070366047488108;
w3=0.107159220467172;
w4=0.139570677926154;
w5=0.166269205816994;
w6=0.186161000115562;
w7=0.198431485327111;
w8=0.202578241925561;
    
R11FP=zeros(2*N+1,2*N+1,L);
R12FP=zeros(2*N+1,2*N+1,L);
R21FP=zeros(2*N+1,2*N+1,L);
R22FP=zeros(2*N+1,2*N+1,L);
r21FP=zeros(2*N+1,2*N+1,L);
ZFP=zeros(2*N+1,2*N+1,L);
epsmn=zeros(4*N+1,L);
 
        
eps=zeros(4*N+1,15*nlayer);
Ieps=zeros(4*N+1,15*nlayer);
Mtm11=zeros(2*N+1,2*N+1,15*nlayer);
Mtm12=zeros(2*N+1,2*N+1,15*nlayer);
Mtm21=zeros(2*N+1,2*N+1,15*nlayer);
Mtm22=zeros(2*N+1,2*N+1,15*nlayer);
 
 
P=leg_gen15(M-1,nlayer);% Legendre polynomials
 
Kx=diag(kx);
d=d/L;
I=eye(2*N+1);
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
 
for lcounter=1:L
    k=(lcounter-1)*15*nlayer+1;
    for kk=1:15*nlayer
%         ----------------------------------------------------
%         Computing Fourier expansion of permittivity function
%         ----------------------------------------------------
        eps(:,kk)=(fastftptrain(xt(:,k+kk-1),epst(:,k+kk-1),2*N)).';
        Ieps(:,kk)=(fastftptrain(xt(:,k+kk-1),1./epst(:,k+kk-1),2*N)).';
%         ----------------------------------------------------
        sizext = size(xt);
        
% --------------------------------------------------
% Evaluating [[c^2]] (here C2) and [[c*s]] (here CS)
% --------------------------------------------------
        mm1=1;
        mm2=1;
        if xt(mm1,k+kk-1)==0
            mm1=mm1+1;
        end
        c2t=0;
        cst=0;
        myxt=0.5;
        while mm1<sizext(1)
            if (xt(mm1,k+kk-1)==xt(mm1+1,k+kk-1))
                mm1=mm1+2;
            else
                if (k+kk-2~=0)&&(k+kk<=length(z))&&(xt(mm1,k+kk)-xt(mm1,k+kk-2)~=0)
                    thetaf = -atan((z(k+kk)-z(k+kk-2))/(xt(mm1,k+kk)-xt(mm1,k+kk-2)));
                    c2t(mm2)=cos(thetaf)^2;
                    cst(mm2)=cos(thetaf)*sin(thetaf);
                    myxt(mm2)=xt(mm1,k+kk-1);
                    mm2 = mm2 + 1;
                end
                mm1 = mm1 + 1;
            end
        end
 
 
        if xt(sizext(1),k+kk-1)~=1
            if mm1>sizext(1)
                mm1 = sizext(1);
            end
            if (k+kk-2~=0)&&(k+kk<=length(z))&&(xt(mm1,k+kk)-xt(mm1,k+kk-2)~=0)
                thetaf = -atan((z(k+kk)-z(k+kk-2))/(xt(sizext(1),k+kk)-xt(sizext(1),k+kk-2)));
                c2t(mm2)=cos(thetaf)^2;
                cst(mm2)=cos(thetaf)*sin(thetaf);
                myxt(mm2)=xt(mm1,k+kk-1);
            end
        end
 
        c2coeff=interpexp(c2t,myxt,1);
        cscoeff=interpexp(cst,myxt,1);
            
        
        c2=zeros(1,4*N+1);
        cs=zeros(1,4*N+1);
        nn=length(c2coeff);
        c2(2*N+1:2*N+nn)=c2coeff;
        cs(2*N+1:2*N+nn)=cscoeff;
% --------------------------------------------------------------------
% Following lines are for structures with sharp edges like triangular
% grating
        if max(abs(cscoeff))>5
             myxt2=ones(1,length(myxt)-1);
            for kk1 = 1:length(myxt)-1
                myxt2(kk1) = (myxt(kk1)+myxt(kk1+1))/2;
            end
            cs=fastftptrain(myxt2,cst,2*N);
        end
% --------------------------------------------------------------------
 
        C2=toeplitz(c2(2*N+1:4*N+1),c2(2*N+1:-1:1));
        CS=toeplitz(cs(2*N+1:4*N+1),cs(2*N+1:-1:1));
        
% ---------------------------------------------------------------
 
% --------------------------------------------
% Calculation of Mtm matrix: See Eqs. (24)-(32)
% --------------------------------------------
        eps_M=toeplitz(eps(2*N+1:4*N+1,kk),eps(2*N+1:-1:1,kk));
        Ieps_M=toeplitz(Ieps(2*N+1:4*N+1,kk),Ieps(2*N+1:-1:1,kk));
        I_Ieps=inv(Ieps_M);
        Delta=(eps_M-I_Ieps);
        A=Delta*C2;
        B=Delta*CS;
        G=-i*inv(eps_M-A);
        KxG=Kx*G;
        BG=B*G;
        Mtm11(:,:,kk)=-KxG*B;
        Mtm12(:,:,kk)=(KxG*Kx/k0+i*k0*I);
        Mtm21(:,:,kk)=k0*(i*A+i*I_Ieps+BG*B);
        Mtm22(:,:,kk)=-BG*Kx;
    end
% --------------------------------------------
 
% --------------------------
% Gauss-Legendre integration
% --------------------------
        for m1=1:M
            for m2=1:m1
                Pm1m2=P(m1,:).*P(m2,:)/nlayer;
                sum1=0;
                sum2=0;
                sum3=0;
                sum4=0;
                for kkk=1:15:15*nlayer-14
                    sum1=sum1+w1*(Pm1m2(kkk)*Mtm11(:,:,kkk)+Pm1m2(kkk+14)*Mtm11(:,:,kkk+14))+w2*(Pm1m2(kkk+1)*Mtm11(:,:,kkk+1)+Pm1m2(kkk+13)*Mtm11(:,:,kkk+13))+w3*(Pm1m2(kkk+2)*Mtm11(:,:,kkk+2)+Pm1m2(kkk+12)*Mtm11(:,:,kkk+12))+w4*(Pm1m2(kkk+3)*Mtm11(:,:,kkk+3)+Pm1m2(kkk+11)*Mtm11(:,:,kkk+11))+w5*(Pm1m2(kkk+4)*Mtm11(:,:,kkk+4)+Pm1m2(kkk+10)*Mtm11(:,:,kkk+10))+w6*(Pm1m2(kkk+5)*Mtm11(:,:,kkk+5)+Pm1m2(kkk+9)*Mtm11(:,:,kkk+9))+w7*(Pm1m2(kkk+6)*Mtm11(:,:,kkk+6)+Pm1m2(kkk+8)*Mtm11(:,:,kkk+8))+w8*(Pm1m2(kkk+7)*Mtm11(:,:,kkk+7));
                    sum2=sum2+w1*(Pm1m2(kkk)*Mtm12(:,:,kkk)+Pm1m2(kkk+14)*Mtm12(:,:,kkk+14))+w2*(Pm1m2(kkk+1)*Mtm12(:,:,kkk+1)+Pm1m2(kkk+13)*Mtm12(:,:,kkk+13))+w3*(Pm1m2(kkk+2)*Mtm12(:,:,kkk+2)+Pm1m2(kkk+12)*Mtm12(:,:,kkk+12))+w4*(Pm1m2(kkk+3)*Mtm12(:,:,kkk+3)+Pm1m2(kkk+11)*Mtm12(:,:,kkk+11))+w5*(Pm1m2(kkk+4)*Mtm12(:,:,kkk+4)+Pm1m2(kkk+10)*Mtm12(:,:,kkk+10))+w6*(Pm1m2(kkk+5)*Mtm12(:,:,kkk+5)+Pm1m2(kkk+9)*Mtm12(:,:,kkk+9))+w7*(Pm1m2(kkk+6)*Mtm12(:,:,kkk+6)+Pm1m2(kkk+8)*Mtm12(:,:,kkk+8))+w8*(Pm1m2(kkk+7)*Mtm12(:,:,kkk+7));
                    sum3=sum3+w1*(Pm1m2(kkk)*Mtm21(:,:,kkk)+Pm1m2(kkk+14)*Mtm21(:,:,kkk+14))+w2*(Pm1m2(kkk+1)*Mtm21(:,:,kkk+1)+Pm1m2(kkk+13)*Mtm21(:,:,kkk+13))+w3*(Pm1m2(kkk+2)*Mtm21(:,:,kkk+2)+Pm1m2(kkk+12)*Mtm21(:,:,kkk+12))+w4*(Pm1m2(kkk+3)*Mtm21(:,:,kkk+3)+Pm1m2(kkk+11)*Mtm21(:,:,kkk+11))+w5*(Pm1m2(kkk+4)*Mtm21(:,:,kkk+4)+Pm1m2(kkk+10)*Mtm21(:,:,kkk+10))+w6*(Pm1m2(kkk+5)*Mtm21(:,:,kkk+5)+Pm1m2(kkk+9)*Mtm21(:,:,kkk+9))+w7*(Pm1m2(kkk+6)*Mtm21(:,:,kkk+6)+Pm1m2(kkk+8)*Mtm21(:,:,kkk+8))+w8*(Pm1m2(kkk+7)*Mtm21(:,:,kkk+7));
                    sum4=sum4+w1*(Pm1m2(kkk)*Mtm22(:,:,kkk)+Pm1m2(kkk+14)*Mtm22(:,:,kkk+14))+w2*(Pm1m2(kkk+1)*Mtm22(:,:,kkk+1)+Pm1m2(kkk+13)*Mtm22(:,:,kkk+13))+w3*(Pm1m2(kkk+2)*Mtm22(:,:,kkk+2)+Pm1m2(kkk+12)*Mtm22(:,:,kkk+12))+w4*(Pm1m2(kkk+3)*Mtm22(:,:,kkk+3)+Pm1m2(kkk+11)*Mtm22(:,:,kkk+11))+w5*(Pm1m2(kkk+4)*Mtm22(:,:,kkk+4)+Pm1m2(kkk+10)*Mtm22(:,:,kkk+10))+w6*(Pm1m2(kkk+5)*Mtm22(:,:,kkk+5)+Pm1m2(kkk+9)*Mtm22(:,:,kkk+9))+w7*(Pm1m2(kkk+6)*Mtm22(:,:,kkk+6)+Pm1m2(kkk+8)*Mtm22(:,:,kkk+8))+w8*(Pm1m2(kkk+7)*Mtm22(:,:,kkk+7));
                end
                Q11(:,:,m1,m2)=sum1;
                Q12(:,:,m1,m2)=sum2;
                Q21(:,:,m1,m2)=sum3;
                Q22(:,:,m1,m2)=sum4;
                Q11(:,:,m2,m1)=Q11(:,:,m1,m2);
                Q12(:,:,m2,m1)=Q12(:,:,m1,m2);
                Q21(:,:,m2,m1)=Q21(:,:,m1,m2);
                Q22(:,:,m2,m1)=Q22(:,:,m1,m2);                                              
            end
        end
% ----------------------------
 
% -----------------------------------------------
% Construction of "Q" matrix: see Eqs. (A5)-(A10)
% -----------------------------------------------
    for m1=1:M-1
        for m2=1:M
%           Eq. (A5):
            Q(m1+(0:(M-1):2*N*(M-1)),m2+(0:M:2*N*M))=Q11(:,:,m1,m2);
%           Eq. (A6):
            Q(m1+(0:(M-1):2*N*(M-1)),m2+(0:M:2*N*M)+(2*N+1)*M)=Q12(:,:,m1,m2);
%           Eq. (A7):
            Q(m1+(0:(M-1):2*N*(M-1))+(2*N+1)*(M-1),m2+(0:M:2*N*M))=Q21(:,:,m1,m2);        
%           Eq. (A8):
            Q(m1+(0:(M-1):2*N*(M-1))+(2*N+1)*(M-1),m2+(0:M:2*N*M)+(2*N+1)*M)=Q22(:,:,m1,m2);
        end
    end
    for m=1:M-1
        for j=-N:N
            for n=m+1+(j+N)*M:2:(j+N)*M+M
%               Eq. (A9):
                Q(m+(j+N)*(M-1),n)=4/d;
%               Eq. (A10):
                Q(m+(j+N)*(M-1)+(2*N+1)*(M-1),n+(2*N+1)*M)=4/d;
            end
        end
    end
% ----------------------------
 
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
        epsmn(:,lcounter)=eps(:,15*nlayer);
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
Uy(N+1,1)=Uy(N+1,1)+k1*cos(alpha)*ux-kx(N+1)*uz;%kz1(N+1)*ux-kx(N+1)*uz;%%
Uy(:,1)=Uy(:,1)/k0;
Sz(:,1) = -Kx*Uy(:,1)/k1;
 
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
 
for j = 2:L+1
    eps_M=toeplitz(epsmn(2*N+1:4*N+1,j-1),epsmn(2*N+1:-1:1,j-1));
    Sz(:,j) = -inv(eps_M)*Kx*Uy(:,j)/k0;
end
 
 
xlen=(L+1)*Nperiod;
x = (0:xlen)*Nperiod/(xlen+1);
z = linspace(0,d*L,L+1);
Hy = zeros(xlen+1,L+1);
Ex = zeros(xlen+1,L+1);
Ez = zeros(xlen+1,L+1);
eta0 = sqrt(4*pi*1e-7/(8.85e-12));
 
for j=-N:N
      Ex=Ex+(exp(-i*kx(j+N+1)*x).')*Sx(j+N+1,:);
      Hy=Hy+(exp(-i*kx(j+N+1)*x).')*Uy(j+N+1,:)/eta0;
      Ez=Ez+(exp(-i*kx(j+N+1)*x).')*Sz(j+N+1,:);
end
