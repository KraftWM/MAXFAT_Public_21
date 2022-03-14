function [Pz_WS_stuetz,Pz_WSD_stuetz,d,f0025] = Pz_Woehlerlinie(Rm,wsgruppe)
% Bestimmt die Stuetzstelle fuer Praj nach FKM Richtlinie Nichtlinear
% Gleichung 2.8-20 & Folgende
% Und rechnet Sie in Pz Stuetzstellen um
% !!! Ohne Sicherheitsfaktor
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Guss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = Hoechstfester Stahl
% Rm       - Zugefestigkeit
% -------------------------------------------------------------------------

% ... PWL
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    apz = 1.173;
    bpz = 1;
    apd = 3.33e-5;
    bpd =  1.55;
    d = -0.56;
    f0025 = 0.39;
    fakD = 1.2414;
    fak = 1.0996;

elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    apz = 10.03;
    bpz = 0.695;
    apd = 5.16e-6;
    bpd = 1.63;
    d = -0.66;
    f0025 = 0.4;
    fakD = 1.2414;
    fak = 1.0768;
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    msg = 'Falsche Werstoffgrupppe uebergeben, Rechnung wird abgebrochen';
    error(msg)  
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    apz = 101.7;
    bpz = 0.26;
    apd = 5.18e-7;
    bpd = 2.04;
    d = -0.61;
    f0025 = 0.36;
    fakD = 1.2414;
    fak = 0.9574;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'Hoechstfester Stahl')
    apz = 0.85;
    bpz = 0.98;
    apd = 4.25*10^-5;
    bpd = 1.44;
    d = -0.56;
    f0025 = 0.31; 
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Rechnung wird abgebrochen';
    error(msg)    
end 

% ... berechne Stuetzstellen
Pz_WS_stuetz = fak * apz * Rm^bpz;
Pz_WSD_stuetz = fakD * apd * Rm^bpd;

end