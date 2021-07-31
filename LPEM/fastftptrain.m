function b=fastftptrain(xt,epst,M)

% fastftptrain(xt,epst,M) calculate 2*M+1 Fourier coefficients of a 
% pulse-like function that its value between xt(i-1) and xt(i) is epst(i).
% The following figure shows one period of a pulse-like function. Notice
% that the period must be normalized to 1 ( 0<xt(i)<1 ).
%                                 epst(4)
%            epst(2)            -----------
%           --------           |           |
%          |        |  epst(3) |           |
%  epst(1) |         ----------            | epst(5)
% ---------                                 --------
% 0      xt(1)    xt(2)      xt(3)       xt(4)     1
% 
% Length(epst) should be length(xt)+1.
% Vector 'b' is the output of function and b(n+M+1) is the n'th Fourier
% coefficient.

b=zeros(1,2*M+1);
lxt=length(xt);
if length(epst) <= lxt
    warning('Length of epst should be more than length of xt')
end
if max(xt)>1 || min(xt)<0
    warning('Period must be normalized to 1 and 0 < xt(i) < 1')
end
for n=-M:M
    if n~=0
        x1=0;
        x2=xt(1);        
        b(n+M+1)=i*epst(1)*(exp(-i*2*pi*n*x2)-exp(-i*2*pi*n*x1))/(2*pi*n);
        for j=2:lxt
            x1=xt(j-1);
            x2=xt(j);
            b(n+M+1)=b(n+M+1)+i*epst(j)*(exp(-i*2*pi*n*x2)-exp(-i*2*pi*n*x1))/(2*pi*n);
        end
        j = lxt;
        x1=xt(j);
        x2=1;        
        b(n+M+1)=b(n+M+1)+i*epst(j+1)*(exp(-i*2*pi*n*x2)-exp(-i*2*pi*n*x1))/(2*pi*n);        
    end
end
x1=0;
x2=xt(1);
b(M+1)=(x2-x1)*epst(1);
for j=2:lxt
    x1=xt(j-1);
    x2=xt(j);
    b(M+1)=b(M+1)+(x2-x1)*epst(j);
end
j = lxt;
x1=xt(j);
x2=1; 
b(M+1)=b(M+1)+(x2-x1)*epst(j+1);