function a=interpexp(f,x,T)
% a=interpexp(f,x,T) --->
% f(x)=a(1)+a(2)*exp(i*2*pi*x/T)+...+a(n)*exp(i*2*pi*(n-1)*x/T)
% f,x are row vectors of the same length
n=length(x);
A=zeros(n);
for j=1:n
    A(:,j)=exp(i*2*pi*(j-1)*x/T).';
end
a=A\(f.');