function [Pfs_WS_stuetz,Pfs_WSD_stuetz,d1,d2,kfs,f0025] = Pfs_WoehlerlinieV3(Rm,wsgruppe)
% Bestimmt die Stuetzstelle fuer Pfs, abschätzung Wächter 
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
    apz = 5.16*10^(-4);
    bpz = 0.5;
    apd = 5.33e-05;
    bpd =  0.5846;
    d1 = -0.45;
    d2 = -0.27;
    f0025 = 0.58;
    ak = 0.4625;
    bk = -0.897;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    apz = 3.88*10^(-4);
    bpz = 0.5;
    apd = 5.33e-05;
    bpd =  0.5846;
    d1 = -0.47;
    d2 = -0.28;
    f0025 = 0.61;
    ak = 0.9006;
    bk = -0.982;
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    apz = 5.7*10^(-4);
    bpz = 0.52;
    apd = 5.33e-05;
    bpd =  0.5846;
    d1 = -0.21;
    d2 = -0.36;
    f0025 = 0.52;
    ak = 0.1823;
    bk = -0.742;
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Werte für Stahl werden verwendet';
    warning(msg)  
    apz = 5.4*10^(-4);
    bpz = 0.5;
    apd = 5.33e-05;
    bpd =  0.5846;
    d1 = -0.35;
    d2 = -0.31;
    f0025 = 0.55;
    ak = 0.4625;
    bk = -0.897;
end 

% ... berechne Stuetzstellen
kfs = ak * Rm ^bk;
Pfs_WS_stuetz = apz * Rm^bpz;
Pfs_WSD_stuetz = apd * Rm^bpd;

end