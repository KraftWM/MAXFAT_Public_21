function fwt = SchubwechselfestigkeitFKM(wsgruppe)
% Bestimmt die Schubwechselfestigkeit nach FKM Richtlinie
%
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Stahlguss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = GJS
%                          6 = GJM
%                          7 = GJL
%                          8 = Hoechstfester Stahl
% -------------------------------------------------------------------------
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    fwt = 1/sqrt(3);
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    fwt = 1/sqrt(3);
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    fwt = 1/sqrt(3);
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    fwt = 0.75;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'GJS')
    fwt = 0.65;
elseif any(wsgruppe == 6) || strcmp(wsgruppe,'GJM')
    fwt = 0.75;
elseif any(wsgruppe == 7) || strcmp(wsgruppe,'GJL')
    fwt = 1;
elseif any(wsgruppe == 8) || strcmp(wsgruppe,'Hoechstfester Stahl')
    fwt = 1/sqrt(3);
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Rechnung wird abgebrochen';
    error(msg)
    
end 

end