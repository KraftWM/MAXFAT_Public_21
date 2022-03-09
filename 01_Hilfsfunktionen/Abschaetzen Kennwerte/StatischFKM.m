function [E,nu] = StatischFKM(wsgruppe)
% Bestimmt die statischen Parameter nach FKM Richtlinie Nichtlinear
% Tabelle 2.30
%
% Querdehnzahl geschätzt 
%
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Guss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = Höchstfester Stahl
% -------------------------------------------------------------------------
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    E = 206000;
    nu = 0.3;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    E = 206000;
    nu = 0.3;
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    E = 70000;
    nu = 0.3;
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    E = 70000;
    nu = 0.3;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'Höchstfester Stahl')
    E = 206000;
    nu = 0.3;
else
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)
    
end 

end