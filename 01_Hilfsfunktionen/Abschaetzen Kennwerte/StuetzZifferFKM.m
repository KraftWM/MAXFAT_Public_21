function [n,nst,nbm] = StuetzZifferFKM(Asig,G,Rm,wsgruppe)
% Bestimmt die Sützwirkung nach FKM Richtlinie Nichtlinear
% Gleichung 2.5-27 & Folgende
% Asig     - Hochbeanspruchte Oberfläche
% G        - bezogener Spannungsgradient
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Guss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = Höchstfester Stahl
% Rm       - Zugefestigkeit
% -------------------------------------------------------------------------

% ... statistische Stüzziffer
% Weibull - Exponent
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    kst = 30;
    Rmbm = 680;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    kst = 15;
    Rmbm = 680;
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    kst = 20;
    Rmbm = 270;
elseif any(wsgruppe == 4) || strcmp(wsgruppe,'Aluguss')
    kst = 20;
    Rmbm = 270;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'Höchstfester Stahl')
    kst = 30;
    Rmbm = 680;
else
    msg = 'Falsche Werstoffgrupppe übergeben, Rechnung wird abgebrochen';
    error(msg)    
end 
% Referenzfläche
Aref = 500; %mm^2
% Stützziffer
nst = (Aref/Asig)^(1/kst);

% ... bruchmechanische Stützziffer
kline = 5 * nst + Rm/Rmbm * sqrt( (7.5 + sqrt(G))/(1+0.2*sqrt(G)) );
nbm = (5 + sqrt(G))/kline;
nbm = max([nbm 1]);

% ... Stuetzziffer
n = nst * nbm;


end