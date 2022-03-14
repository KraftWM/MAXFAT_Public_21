function eps = transformstrain(varargin)
% Funktion transformiert Dehnungszustand aus dem Koordinantensystem der
% Kerbe in das lokale Koordinantensystem der Schnittebene. Dabei entsteht 
% im Allgemeinen aus einem ebenen - ein raeumlicher Spannungszustand.
%
% !!! Nebendiagonalelemente werden als Gleitungen behandelt
%
% INPUT:
%    varargin  - Variabler Input
%
% -------------------------------------------------------------------------
% Modus 1: nargin == 3
%    EPS -> Dehnungen im ESZ im Kerbkoordinantensystems
%    phi -> Drehwinkel des lokalen Koord. Winkel zwischen Y,y 
%           angegeben in bogenmaß
%    psi -> Drehwinkel des lokalen Koord. Winkel zwischen Z,z 
%           angegeben im Bogenmaß
%
% OUTPUT:
%   eps -> Dehnungstensor
%
% -------------------------------------------------------------------------
% Modus 2: nargin == 2
%    phi -> Drehwinkel des lokalen Koord. Winkel zwischen Y,y 
%           angegeben im Bogenmaß
%    psi -> Drehwinkel des lokalen Koord. Winkel zwischen Z,z 
%           angegeben im Bogenmaß
%
% OUTPUT:
%   R    -> Rotationstensor fuer Dehnungen (6x4)
% _________________________________________________________________________

% unterscheide variablen input
if nargin == 2
    phi = varargin{1,1};
    psi = varargin{1,2};
elseif nargin == 3
    EPS = varargin{1,1};
    phi = varargin{1,2};
    psi = varargin{1,3};
end

% Rotationsmatrix
R11 = cos(phi) * cos(psi);
R12 = sin(phi) * cos(psi);
R13 = sin(psi);
R21 = -sin(phi);
R22 = cos(phi);
R23 = 0;
R31 = -cos(phi) * sin(psi);
R32 = -sin(phi) * sin(psi);
R33 = cos(psi);

R = [...
         R11^2 ,    R12^2 ,    R13^2 ,     R11*R12      ;...
         R21^2 ,    R22^2 ,    R23^2 ,     R21*R22      ;...
         R31^2 ,    R32^2 ,    R33^2 ,     R31*R32      ;...
     2*R11*R21 , 2*R12*R22, 2*R13*R23, R11*R22 + R12*R21;...
     2*R21*R31 , 2*R22*R32, 2*R23*R33, R21*R32 + R22*R31;...
     2*R31*R11 , 2*R32*R12, 2*R33*R13, R31*R12 + R32*R11;...
    ];


% Transformation der Dehnungen
if nargin == 3
    eps = R * EPS;
elseif nargin == 2
    eps = R;
end

end % Ende Funktion