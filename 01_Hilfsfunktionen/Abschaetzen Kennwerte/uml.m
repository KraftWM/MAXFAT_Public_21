function [sf, ef, tf, gf, b, c, ND, Ks, ns] = uml(E, Rm, grp)
   
% Funktion berechnet Materialparameter nach dem Unified Material Law
% aus [1] 'Abschätzung von zyklischen Werkstoffkennwerten' (A. Hatscher)
% 
% INPUT:
% E         = Elastizitätsmodul
% Rm        = mittlere Zugfestigkeit
% grp       = Werkstoffgruppe nach [1], Tab.1, UML
%               grp = 1 : un- und niedriglegierte Stähle
%               grp = 3 : Al- und Ti-Legierungen     (gruppen nach
%                                                      RiliNili)
% OUTPUT:
% sf        = zyklischer Schwingfestigkeitskoeffizient
% ef        = zyklischer Duktilitätskoeffizient
% tf        = zyklischer Schubfestigkeitskoeffizient
% gf        = zyklischer Gleitungskoeffizient
% b         = zyklischer Schwingfestigkeitsexponent
% c         = zyklischer Duktilitätsexponent
% ND        = Abknickpunkt Dauerfestigkeit
% Ks,ns     = Parameter Ramberg Osgood
%__________________________________________________________________________

% ... Check Gültigkeit des inputs
if Rm < 110
    Rm = 110;
elseif Rm > 2300
    Rm = 2300;
end

% ... bestimme psi
if Rm/E <= 0.003
    psi = 1;
else
    psi = 1.375 - 125*Rm/E;
end

% ... Unterscheide WS Gruppe
if grp == 1 || grp == 2 || grp == 4
    ND = 1e6;%5e5;
    b = -0.087; c = -0.58;
    sf = 1.5*Rm;
    ef = 0.59*psi;
    Ks = 1.65*Rm;
    ns = 0.15;
elseif grp == 3
    ND = 1e6;
    b = -0.095; c = -0.69;
    sf = 1.67*Rm;
    ef = 0.35;
    Ks = 1.61 * Rm;
    ns = 0.11;
end

% ... Gleitungswöhlerlinie abschätzen
tf = sf/(sqrt(3));
gf = ef*(sqrt(3));


end