function KPR = RauheitseinflussFKM(wsgruppe,Rz,Rm)
% Rauheitseinfluss nach FKM Richtlinie nichtlinear
% Tabelle 2.35 & Gleichung 2.8-37
%
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Guss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = Hoechstfester Stahl
% Rz       - Rauheit in Mikrometer 
% Rm       - Zugefestigkeit in MPa
% -------------------------------------------------------------------------

% Abfangen glatte Oberflaeche
if Rz == 0
    KPR = 1;
    return;
end


if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    arp = 0.27;
    brp = 0.43;
    RmNmin = 400;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    arp = 0.25;
    brp = 0.42;
    RmNmin = 400;
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    msg = 'Falsche Werstoffgrupppe uebergeben, Rechnung wird abgebrochen';
    error(msg)  
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    arp = 0.27;
    brp = 0.43;
    RmNmin = 133;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'Hoechstfester Stahl')
    arp = 0.27;
    brp = 0.43;
    RmNmin = 400; 
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Rechnung wird abgebrochen';
    error(msg)    
end 

KPR = (1-arp*log10(Rz)*log10(2*Rm/RmNmin))^brp;

end