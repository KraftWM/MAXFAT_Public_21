function [Praj_WS_stuetz,Praj_WSD_stuetz,d,f0025] = Praj_Woehlerlinie_FKM(Rm,wsgruppe)
% Bestimmt die Stuetzstelle fuer Praj nach FKM Richtlinie Nichtlinear
% Gleichung 2.8-20 & Folgende 
% !!! Ohne Sicherheitsfaktor
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

% ... PWL
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    apz = 1.173;
    bpz = 1;
    apd = 3.33e-5;
    bpd =  1.55;
    d = -0.56;
    f0025 = 0.35;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    apz = 230;
    bpz = 0.35;
    apd = 5.16e-6;
    bpd = 1.63;
    d = -0.72;
    f0025 = 0.34;
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    apz = 213;
    bpz = 243;
    apd = 5.18e-7;
    bpd = 2.04;
    d = -0.61;
    f0025 = 0.33;
elseif any(wsgruppe == 8) || strcmp(wsgruppe,'Hoechstfester Stahl')
    apz = 0.85;
    bpz = 0.98;
    apd = 4.25*10^-5;
    bpd = 1.44;
    d = -0.56;
    f0025 = 0.31; 
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Werte für Stahl übernommen';
    warning(msg)
    apz = 10;
    bpz = 0.826;
    apd = 3.33e-5;
    bpd =  1.55;
    d = -0.63;
    f0025 = 0.35;
end 

% ... berechne Stuetzstellen
Praj_WS_stuetz = apz * Rm^bpz;
Praj_WSD_stuetz = apd * Rm^bpd;

end