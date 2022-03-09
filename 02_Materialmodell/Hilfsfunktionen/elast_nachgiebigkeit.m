function D = elast_nachgiebigkeit(E,nu,ntens,ndi)
% Funktion berechnet Steifigkeitstensor, je nach Spannungsszustand
if ntens == 6 % 3D
    D = 1/E .*[1,   -nu,   -nu,          0,          0, 0;...
               -nu,   1,   -nu,          0,          0, 0;...
               -nu, -nu,     1,          0,          0, 0;
                 0,    0,    0,        2*(1+nu),     0, 0;
                 0,    0,    0,          0,   2*(1+nu), 0;...
                 0,    0,    0,          0,          0,2*(1+nu)];
elseif ntens == 3 && ndi == 2 % ESZ
    D = 1/E .*[1,   -nu, 0;...
               -nu,   1, 0;...
                 0,   0, 2*(1+nu)];
elseif ntens == 2 && ndi == 1 % sigma-tau
    D = [1/E,0;...
         0,2*(1+nu)/E];
elseif ntens == 1 % 1D
    D = 1/E;
else % Damit Compiler in Coder nicht rumzickt
    D = 1;
end