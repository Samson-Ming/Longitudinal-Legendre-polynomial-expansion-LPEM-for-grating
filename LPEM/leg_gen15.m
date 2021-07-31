function P=leg_gen15(M,n)
% P=leg_gen15(M,n) generates Legendre functions. Intervals are chosen to be 
% compatible with Gauss-Legendre(15 points) integration algorithm.
% The output P is a (M,15*n) matrix and P(m+1,:) is the m'th Legendre
% function.
% n: number of sections. This implies that we have used 15*n points.
% M: number of Legendre polynomials

zlen=n*15;
delta=2/n;
P=zeros(M,zlen);
a=-1;
b=1;
a1=a-delta;
x1=-0.987992518020485;
x2=-0.937273392400706;
x3=-0.848206583410427;
x4=-0.724417731360170;
x5=-0.570972172608539;
x6=-0.394151347077563;
x7=-0.201194093997435;
x8=0;
x9=-x7;
x10=-x6;
x11=-x5;
x12=-x4;
x13=-x3;
x14=-x2;
x15=-x1;
for j=1:n
    a2=a1+delta;
    b2=a2+delta;
    a1=a2;
    x(1+15*(j-1))=((b2-a2)*x1+(b2+a2))/2;
    x(2+15*(j-1))=((b2-a2)*x2+(b2+a2))/2;
    x(3+15*(j-1))=((b2-a2)*x3+(b2+a2))/2;
    x(4+15*(j-1))=((b2-a2)*x4+(b2+a2))/2;
    x(5+15*(j-1))=((b2-a2)*x5+(b2+a2))/2;
    x(6+15*(j-1))=((b2-a2)*x6+(b2+a2))/2;
    x(7+15*(j-1))=((b2-a2)*x7+(b2+a2))/2;
    x(8+15*(j-1))=((b2-a2)*x8+(b2+a2))/2;
    x(9+15*(j-1))=((b2-a2)*x9+(b2+a2))/2;
    x(10+15*(j-1))=((b2-a2)*x10+(b2+a2))/2;
    x(11+15*(j-1))=((b2-a2)*x11+(b2+a2))/2;
    x(12+15*(j-1))=((b2-a2)*x12+(b2+a2))/2;
    x(13+15*(j-1))=((b2-a2)*x13+(b2+a2))/2;
    x(14+15*(j-1))=((b2-a2)*x14+(b2+a2))/2;
    x(15+15*(j-1))=((b2-a2)*x15+(b2+a2))/2;
end
P(1,:)=1;
P(2,:)=x;
for n=2:M
    P(n+1,:)=2*x.*P(n,:)-P(n-1,:)-(x.*P(n,:)-P(n-1,:))/n;
end
