function [tx, ty, tz] = planeprojection(tensor,ex,ey,ez)
%
% Funktion Projeziert gegebenen 3D Tensorzustand in eine Ebene, die 
% ex als Normalenvektor besitzt
%
% INPUT:
% tensor   -> 3D Tensor(verlauf) in Voigt notation, gedacht für 
%             Spannungs- & Dehnungstensoren. Für Dehnungstensoren doppelte 
%             Nebendiagonalelemente
% ex       -> Normaleneinheitsvektore der Ebene
% ey       -> Einheitsvektor in der Ebene
% ez       -> Einheitsvektor in der Ebene (ex,ey,ez) als orthonormalsystem
%
% OUTPUT:
% tx       -> komponente in x-richtung (Bepsiel sigxx, epsxx)
% ty       -> komponente in y-richtung
% tz       -> komponente in z-richtung
%
%
% Notationen
%         sigxx              epsxx
%         sigyy              epsyy
%  sig =  sigzz       eps =  epszz
%         sigxy             2epsxy
%         sigyz             2epsyz
%         sigxz             2epsxz
%
% !!!! Kommentare für Spannungen !!!!!
%__________________________________________________________________________

% Schnittspannung nach Chauchy Theorem (t_i = sig_ij n_j)

t = [...
    tensor(1,:) * ex(1) + tensor(4,:) * ex(2) + tensor(6,:) * ex(3); ...
    tensor(4,:) * ex(1) + tensor(2,:) * ex(2) + tensor(5,:) * ex(3); ...
    tensor(6,:) * ex(1) + tensor(5,:) * ex(2) + tensor(3,:) * ex(3); ...
    ]; 

% Normalen Komponente
tx = ex' * t;
% y Komponente
ty = ey' * t;
% z Komponente
tz = ez' * t;

end % Ende Funktion