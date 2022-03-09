function C = elast_steifigkeit(E,nu,ntens,ndi)
% Funktion berechnet Steifigkeitstensor, je nach Spannungsszustand
if ntens == 6 % 3D
    C = E/((1+nu)*(1-2*nu)).*[1-nu,   nu,   nu,          0,          0, 0;...
                            nu, 1-nu,   nu,          0,          0, 0;...
                            nu,   nu, 1-nu,          0,          0, 0;
                             0,    0,    0, (1-2*nu)/2,          0, 0;
                             0,    0,    0,          0, (1-2*nu)/2, 0;...
                             0,    0,    0,          0,           0,(1-2*nu)/2];
elseif ntens == 3 && ndi == 2 % ESZ
    C = E/(1- nu^2) * [1 ,nu,0; nu,1 ,0;0 ,0 ,(1-nu)/2];
elseif ntens == 2 && ndi == 1 % sigma-tau
    C = [E,0;...
         0,E/(2*(1+nu))];
elseif ntens == 1 % 1D
    C = E;
else % Dummy für Coder
    C = 1;
end