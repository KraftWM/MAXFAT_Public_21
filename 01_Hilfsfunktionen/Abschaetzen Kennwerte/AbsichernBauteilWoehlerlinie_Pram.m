function fram = AbsichernBauteilWoehlerlinie_Pram(PA,n,KRP,f0025)
% Bestimmt Faktor zum Absichern der Bauteilwöhlerlinie von PRAM
% nach Richtlinie Nichtlinear
%
% INPUT:
% PA    - Ausfallwahrscheinlichkeit
% n     - Stützziffer
% KPR   - Rauheitsfaktor
% f0025 - Faktor zur Absicherung der Werkstoff-WL bei
%         Ausfallwahrscheinlichkeit PA = 2.5%
% -------------------------------------------------------------------------


% Absichern der Bauteilwöhlerlinie nach FKM-Richtlinie Nichtlinear für Pram
% Tabelle 2.36
if abs(PA - 2.3 * 10^(-1)) < 1e-5
    gamma_M=1.1;
elseif abs(PA - 10^(-3)) < 1e-5
    gamma_M=1.1;
elseif abs(PA - 7.2 * 10^(-5)) < 1e-5
    gamma_M=1.2;
elseif abs(PA - 10^(-5)) < 1e-5
    gamma_M=1.3;
elseif abs(PA - 0.5) < 1e-5 
    gamma_M = 1;
else
    msg = 'berechnen ohne Absicherung der PRAM Wöhlerlinie';
    warning(msg);
    gamma_M = 1;
end

% Zusammenfassen aller Faktoren
if abs(PA - 0.5) < 1e-5
    fram = gamma_M/(n*KRP);
else
    fram = gamma_M/(f0025*n*KRP);
end