function [Kstrich,nstrich] = RambergOsgoodFKM(wsgruppe,Rm)
% Bestimmt den zyklische Parameter nach FKM Richtlinie Nichtlinear
% Gleichung 2.6-10 & Tabelle 2.17
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Stahlguss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = GJS
%                          6 = GJM
%                          7 = GJL
%                          8 = Hoechstfester Stahl
% Rm       - Zugefestigkeit
% -------------------------------------------------------------------------
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    nstrich = 0.187;
    asig = 3.1148;
    aeps = 1033;
    bsig = 0.897;
    beps = -1.235;
    epsgrenz = 0.338;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    nstrich = 0.176;
    asig = 1.732;
    aeps = 0.847;
    bsig = 0.982;
    beps = -0.181;
    epsgrenz = 1e40;
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    nstrich = 0.128;
    asig = 9.12;
    aeps = 895.9;
    bsig = 0.742;
    beps = -1.183;
    epsgrenz = 1e40;
elseif any(wsgruppe == 8) || strcmp(wsgruppe,'Hoechstfester Stahl')
    asig = 2.66;
    aeps = 1400;
    bsig = 0.895;
    beps = -1.235;
    epsgrenz = 0.099;
    nstrich = 0.085;
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Werte für Stahl werden verwendet';
    warning(msg)
    nstrich = 0.187;
    asig = 3.1148;
    aeps = 1033;
    bsig = 0.897;
    beps = -1.235;
    epsgrenz = 0.338;
    
end 

Zahler = asig * Rm^bsig;
Nenner = (min([epsgrenz aeps*Rm^beps]))^nstrich;
Kstrich = Zahler/Nenner;
end