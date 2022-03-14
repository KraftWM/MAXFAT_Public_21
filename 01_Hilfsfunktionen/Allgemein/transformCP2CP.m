function R = transformCP2CP(phi1,phi2,psi1,psi2,sigeps_flag)
% Funktion gibt Rotationsmatrix fuer Tensoren in Voigt Notation zurueck
% Rotiert wird zwischen zwei kritischen Ebenen.
%
% K0   -> Kerb Koordinatensysten (Ebener Spannungszustand)
% K1   -> Koordinatensystem krtische Ebene 1 (3D Spannungszustand), bekannt
% K2   -> Koordinatensystem krtische Ebene 2 (3D Spannungszustand)
%
% Rotationsmatrix rotiert von K1 -> K2
%  2_SIG = R * 1_SIG
%
% Voigt Notation fuer Spannungen und Dehnungen in der kritischen Ebene:
%         sigxx              epsxx
%         sigyy              epsyy
%  SIG =  sigzz       EPS =  epszz
%         sigxy             2epsxy
%         sigyz             2epsyz
%         sigxz             2epsxz
%
% INPUT:
%    phi1, psi1   - Rotationswinkel von K0->K1 im Bogenmaß
%    phi2, psi2   - Rotationswinkel von K0->K2 im Bogenmaß  
%    sigeps_flag  - Unterscheide Rotation von Spannungen und Dehnngen
%                   (wegen verdopplung der ND Elemente)
%                   1 - Dehnungsrotation
%                   0 - Spannungsrotation
% OUTPUT:
%   R    -> Rotationstensor (6x6)
%
% -------------------------------------------------------------------------

% ... Winkelinkrement
dphi = phi2 - phi1;

% ... Komponenten Rotaionstensor
R11 = cos(psi1)*cos(psi2)*cos(dphi)+sin(psi1)*sin(psi2);
R12 = cos(psi2)*sin(dphi);
R13 = -cos(psi2)*sin(psi1)*cos(dphi)+sin(psi2)*cos(psi1);

R21 = -cos(psi1)*sin(dphi);
R22 = cos(dphi);
R23 = sin(psi1)*sin(dphi);

R31 = -cos(psi1)*sin(psi2)*cos(dphi) + cos(psi2)*sin(psi1);
R32 = -sin(psi2)*sin(dphi);
R33 = sin(psi2)*sin(psi1)*cos(dphi)+cos(psi2)*cos(psi1);

% ... Rotationsmatrix fuer Voigt Notation
% Dehnungstransformation
if sigeps_flag
    R = [    R11^2      R12^2      R13^2          R11*R12          R12*R13          R11*R13;...
             R21^2      R22^2      R23^2          R22*R21          R22*R23          R21*R23;...
             R31^2      R32^2      R33^2          R31*R32          R33*R32          R33*R31;...
         2*R21*R11  2*R22*R12  2*R23*R13  R21*R12+R22*R11  R22*R13+R23*R12  R21*R13+R23*R11;...
         2*R21*R31  2*R22*R32  2*R23*R33  R31*R22+R32*R21  R32*R23+R33*R22  R31*R23+R33*R21;...
         2*R31*R11  2*R32*R12  2*R33*R13  R31*R12+R32*R11  R32*R13+R33*R12  R31*R13+R33*R11];
% Spannungstransformation
else
    R = [  R11^2    R12^2    R13^2        2*R11*R12        2*R12*R13        2*R11*R13;...
           R21^2    R22^2    R23^2        2*R22*R21        2*R22*R23        2*R21*R23;...
           R31^2    R32^2    R33^2        2*R31*R32        2*R33*R32        2*R33*R31;...
         R21*R11  R22*R12  R23*R13  R21*R12+R22*R11  R22*R13+R23*R12  R21*R13+R23*R11;...
         R21*R31  R22*R32  R23*R33  R31*R22+R32*R21  R32*R23+R33*R22  R31*R23+R33*R21;...
         R31*R11  R32*R12  R33*R13  R31*R12+R32*R11  R32*R13+R33*R12  R31*R13+R33*R11];
    
end
end