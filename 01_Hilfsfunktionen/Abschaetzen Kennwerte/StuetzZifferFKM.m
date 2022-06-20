function [n,nst,nbm] = StuetzZifferFKM(Asig,G,Rm,wsgruppe)
% Bestimmt die Suetzwirkung nach FKM Richtlinie Nichtlinear
% Gleichung 2.5-27 & Folgende
% Asig     - Hochbeanspruchte Oberflaeche
% G        - bezogener Spannungsgradient
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

% ... statistische Stuezziffer
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
    kst = 10;
    Rmbm = Inf;
elseif any(wsgruppe == 5) || strcmp(wsgruppe,'GJS')
    kst = 10;
    Rmbm = Inf;
elseif any(wsgruppe == 6) || strcmp(wsgruppe,'GJM')
    kst = 10;
    Rmbm = Inf;
elseif any(wsgruppe == 7) || strcmp(wsgruppe,'GJL')
    kst = 10;
    Rmbm = Inf;
elseif any(wsgruppe == 8) || strcmp(wsgruppe,'Hoechstfester Stahl')
    kst = 30;
    Rmbm = 680;
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Rechnung wird abgebrochen';
    error(msg)    
end 
% Referenzflaeche
Aref = 500; %mm^2
% Stuetzziffer
nst = (Aref/Asig)^(1/kst);

% ... bruchmechanische Stuetzziffer
kline = 5 * nst + Rm/Rmbm * sqrt( (7.5 + sqrt(G))/(1+0.2*sqrt(G)) );
nbm = (5 + sqrt(G))/kline;
nbm = max([nbm 1]);
if wsgruppe >= 4
    nbm = 1;
end

% ... Stuetzziffer
n = nst * nbm;


end