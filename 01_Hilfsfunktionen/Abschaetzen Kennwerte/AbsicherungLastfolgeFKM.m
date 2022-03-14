function gamma_L = AbsicherungLastfolgeFKM(verfahren,sL,Lmax,PL,PA)
% Statistische Absicherung der Lastfolge nach FKM Richtlinie Nichtlinear
% siehe Abschnitt 2.3.2
%
% INPUT:
% verfahren  : 1 - Annahme einer normalverteilten Streuung mit 
%                  Standardabweichung sL
%
%              2 - Annahme einer logarithmisch-normalverteilten Streuung 
%                  mit Standardabweichung LSDs nach FKM-Richtlinie
%
%              3 - Als Pauschalwert bei unbekannter Streuung fuer PL = 2,5%
%
%  sl        - Streuung normal oder log normal Verteilung je nach Verfahren
%  Lmax      - maximalwert Lastfolge;
%  PL        - Auftretenswahrscheinlichkeit der Lastfolge
%  PA        - Ausfallwahrscheinlichkeit 
%
% -------------------------------------------------------------------------

% Berechne beta
if abs(PA - 2.3 * 10^(-1)) < 1e-5
    beta=0.739;
elseif abs(PA - 10^(-3)) < 1e-5
    beta=3.09;
elseif abs(PA - 7.2 * 10^(-5)) < 1e-5
    beta=3.8;
elseif abs(PA - 10^(-5)) < 1e-5
    beta=4.27;
elseif abs(PA - 0.5) < 1e-5 
    gamma_L = 1;
    return;
else
    msg = 'berechnen ohne Lastabsicherung';
    warning(msg);
    gamma_L = 1;
    return;
end

% Unterscheide verfahren
switch verfahren % Normalverteilung
    case 1        
        % Beiwert
        if PL == 0.025
            alpha_L = (0.7*beta-2)*sL;
        elseif PL == 0.5
            alpha_L = 0.7*beta*sL;
        else
            msg = 'falsche Auftretenswahrscheinlichkeit angegeben, es wird mit PL = 50% gerechnet';
            warning(msg)
            alpha_L = 0.7*beta*sL;
        end
        % Sicherheitsfaktor
        gamma_L = (Lmax + alpha_L ) / Lmax;
        gamma_L = max([gamma_L 1]);
    case 2 % Log Normalverteilung
        % Beiwert
        if PL == 0.025
            alpha_L = (0.7*beta-2)*sL;
        elseif PL == 0.5
            alpha_L = 0.7*beta*sL;
        else
            msg = 'falsche Auftretenswahrscheinlichkeit angegeben, es wird mit PL = 50% gerechnet';
            warning(msg)
            alpha_L = 0.7*beta*sL;
        end
        % Sicherheitsfaktor
        gamma_L = 10^alpha_L;
        gamma_L = max([gamma_L 1]);
    case 3 % sonst
        gamma_L = 1.1;
    otherwise
        msg = 'Falsches Verfahren bei Absicherung der Lastfolge angegeben';
        error(msg)
end

end