function [Pfs_WS_stuetz,Pfs_WSD_stuetz,d1,d2] = Pfs_Woehlerlinie(Rm,wsgruppe)
% Bestimmt die Stützstelle für Pfs ... sehr sehr vorläufige Version
% !!! Ohne Sicherheitsfaktor
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Guss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = Höchstfester Stahl
% Rm       - Zugefestigkeit
% -------------------------------------------------------------------------

% ... PWL
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    apz = 1.348;
    bpz = -0.6745;
    apd = 5.33e-05;
    bpd =  0.5846;
    d1 = -0.5043;
    d2 = -0.2506;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)  
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)  
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)  
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'Höchstfester Stahl')
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)  
else
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)    
end 

% ... berechne Stützstellen
Pfs_WS_stuetz = apz * Rm^bpz;
Pfs_WSD_stuetz = apd * Rm^bpd;

end