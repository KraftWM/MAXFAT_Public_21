function ffs = AbsichernBauteilWoehlerlinie_Pfs(PA,n,KRP,f0025)
% Bestimmt Faktor zum Absichern der Bauteilwoehlerlinie von Pfs
% aus Abschlussbericht FKM Projekt "Mehrachsig Örtlich" 
%
% INPUT:
% PA   - Ausfallwahrscheinlichkeit
% n    - Stützziffer
% KPR  - Rauheitsfaktor
% f0025 - Faktor zur Absicherung der Werkstoff-WL bei
%         Ausfallwahrscheinlichkeit PA = 2.5%
% OUTPUT:
% fz  - Sicherheitsfaktor Pfs
% -------------------------------------------------------------------------

% Absichern der Bauteilwoehlerlinie nach FKM-Richtlinie Nichtlinear für Pram
% Tabelle 2.13
if abs(PA - 2.3 * 10^(-1)) < 1e-5
    gamma_M=1.7;
elseif abs(PA - 10^(-3)) < 1e-5
    gamma_M=3.0;
elseif abs(PA - 7.2 * 10^(-5)) < 1e-5
    gamma_M=4.35;
elseif abs(PA - 10^(-5)) < 1e-5
    gamma_M=5.25;
elseif abs(PA - 0.5) < 1e-5 
    gamma_M = 1;
else
    msg = 'berechnen ohne Absicherung der Praj Woehlerlinie';
    warning(msg);
    gamma_M = 1;
end

% Zusammenfassen aller Faktoren
if abs(PA - 0.5) < 1e-5
    ffs = gamma_M/(n*KRP)^2;
else
    ffs = gamma_M/(n*KRP)^2/f0025;
end
