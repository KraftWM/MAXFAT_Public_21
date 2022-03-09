function fraj = AbsichernBauteilWoehlerlinie_Praj(PA,n,KRP,f0025)
% Bestimmt Faktor zum Absichern der Bauteilwöhlerlinie von PRAM
% nach Richtlinie Nichtlinear
%
% INPUT:
% PA   - Ausfallwahrscheinlichkeit
% n    - Stützziffer
% KPR  - Rauheitsfaktor
% f0025 - Faktor zur Absicherung der Werkstoff-WL bei
%         Ausfallwahrscheinlichkeit PA = 2.5%
% -------------------------------------------------------------------------

% Absichern der Bauteilwöhlerlinie nach FKM-Richtlinie Nichtlinear für Pram
% Tabelle 2.13
if abs(PA - 2.3 * 10^(-1)) < 1e-5
    gamma_M=1.2;
elseif abs(PA - 10^(-3)) < 1e-5
    gamma_M=1.2;
elseif abs(PA - 7.2 * 10^(-5)) < 1e-5
    gamma_M=1.45;
elseif abs(PA - 10^(-5)) < 1e-5
    gamma_M=1.7;
elseif abs(PA - 0.5) < 1e-5 
    gamma_M = 1;
else
    msg = 'berechnen ohne Absicherung der PRAM Wöhlerlinie';
    warning(msg);
    gamma_M = 1;
end

% Zusammenfassen aller Faktoren
if abs(PA - 0.5) < 1e-5
    fraj = gamma_M/(n*KRP)^2;
else
    fraj = gamma_M/(n*KRP)^2/f0025;
end
