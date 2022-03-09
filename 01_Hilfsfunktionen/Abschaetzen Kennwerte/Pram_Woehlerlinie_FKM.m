function [Pram_WS_stuetz,Pram_WSD_stuetz,d1,d2,f0025] = Pram_Woehlerlinie_FKM(Rm,wsgruppe)
% Bestimmt die Stützstelle für Pram nach FKM Richtlinie Nichtlinear
% Gleichung 2.5-22 & Folgende 
% !!! Ohne sicherheitsfaktor
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Guss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = Höchstfester Stahl
% Rm       - Zugefestigkeit
% -------------------------------------------------------------------------

% ... PWL
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    apz = 20;
    bpz = 0.587;
    apd = 0.82;
    bpd =  0.92;
    d1 = -0.302;
    d2 = -0.197;
    f0025 = 0.71;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    apz = 25.56;
    bpz = 0.519;
    apd = 0.46;
    bpd = 0.96;
    d1 = -0.289;
    d2 = -0.189;
    f0025 = 0.51;
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)  
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    apz = 16.71;
    bpz = 0.537;
    apd = 0.3;
    bpd = 1;
    d1 = -0.238;
    d2 = -0.167;
    f0025 = 0.61;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'Höchstfester Stahl')
    apz = 18;
    bpz = 0.587;
    apd = 0.73;
    bpd =  0.93;
    d1 = -0.155;
    d2 = -0.145;
    f0025 = 0.65; 
else
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)    
end 

% ... berechne Stützstellen
Pram_WS_stuetz = apz * Rm^bpz;
Pram_WSD_stuetz = apd * Rm^bpd;

end