function a=fastft(y1,M)
% fastft(y1,M) gives 2*M+1 Fourier coefficients of y1.
%In this function it has been assumed that the period of y1 is 1. M is the last
%Fourier coefficient and size(a)=2*M+1 . The m'th Fourier coefficient of y1 
%a(m+M+1)
siz=length(y1);
n=(siz-1)/2;
ef=fft(real(y1))/(2*n+1);
ef2=fft(imag(y1))/(2*n+1);
e = zeros(1, 2*n-1);
e(n : 2*n-1) = ef(1 : n)+i*ef2(1:n);
e(1 : n-1) = conj(ef(n :-1: 2))+i*conj(ef2(n:-1:2));
for m=-M:M
    a(m+M+1)=e(n+m);
end