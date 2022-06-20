function [Pz_WS_stuetz,Pz_WSD_stuetz,d,f0025] = Pz_Woehlerlinie(Rm,wsgruppe)
% Bestimmt die Stuetzstelle fuer Pz
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
    apz = 50;
    bpz = 0.567;
    apd = 3.33e-5;
    bpd =  1.55;
    d = -0.61;
    f0025 = 0.38;

elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    apz = 208;
    bpz = 0.355;
    apd = 5.16e-6;
    bpd = 1.63;
    d = -0.69;
    f0025 = 0.33;

elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    apz = 192;
    bpz = 0.243;
    apd = 5.18e-7;
    bpd = 2.04;
    d = -0.63;
    f0025 = 0.37;

else
    msg = 'Falsche Werstoffgrupppe uebergeben, Werte für Stahl übernommen';
    warning(msg) 
    apz = 50;
    bpz = 0.567;
    apd = 3.33e-5;
    bpd =  1.55;
    d = -0.61;
    f0025 = 0.38;
end 

% ... berechne Stuetzstellen
Pz_WS_stuetz = apz * Rm^bpz;
Pz_WSD_stuetz = apd * Rm^bpd;

end