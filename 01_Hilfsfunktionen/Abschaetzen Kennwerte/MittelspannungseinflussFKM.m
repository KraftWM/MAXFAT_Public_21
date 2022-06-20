function Msig = MittelspannungseinflussFKM(wsgruppe,Rm)
% Bestimmt den Mittelspannungseinfluss nach FKM Richtlinie Nichtlinear
% Gleichung 2.5-61 & Tabelle 2.14
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
    am = 0.35;
    bm = -0.1;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    am = 0.35;
    bm = 0.05;
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    am = 1;
    bm = -0.04;
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    am = 1;
    bm = 0.2;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'GJS')
    am = 0.35;
    bm = 0.08;
elseif any(wsgruppe == 6) || strcmp(wsgruppe,'GJM')
    am = 0.35;
    bm = 0.13;
elseif any(wsgruppe == 7) || strcmp(wsgruppe,'GJL')
    am = 0;
    bm = 0.5;
elseif any(wsgruppe == 8) || strcmp(wsgruppe,'Hoechstfester Stahl')
    am = 0.39;
    bm = -0.36;
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Rechnung wird abgebrochen';
    error(msg)    
end 

Msig = am * Rm * 0.001 + bm;

end