function kappaz = Blochwavenumber(R11,R12,R21,R22,d)
% This function calculates normal component of Bloch wave number by using 
% the R matrix of the structure. Notice that the other component depends on
% the incident angle in which the R matrix has been obtained (kappax =
% 2*pi*n1*sin(alpha)/lambda where n1 is the incident medium's refractive index and
% alpha is incident angle.).
% 
%                       Region 1 (Incident medium)
%                  ------------------                |-----> x
%                   Periodic Region                  |     
%                  ------------------                V
%                       Region 3                     z
% 
% d: Grating thickness (Normalized to the grating period)
% N: Truncation order
% R11, R12, R21 and R22 are sub-matrixes of the R matrix.
% R =[R11 R12
%     R21 R22]
sizeR11 = size(R11);
N = (sizeR11(1) - 1) / 2;
I=eye(2*N+1);
Zero = zeros(2*N+1);
eigT=eig([-I R22;Zero R12],[Zero -R21;I -R11]);
[bikhod Indiceeig ]=min(abs(abs(eigT)-1));
kappaz=abs(angle(eigT(Indiceeig))/d); 
