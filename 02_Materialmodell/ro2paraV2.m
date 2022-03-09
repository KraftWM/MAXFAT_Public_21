function para = ro2paraV2(typ,...
                          E,nu,Kstrich,nstrich,M,...
                          modell,verfahren_flag,...
                          varargin)
% Funktion zum bestimmen der Materialparameter aus zyklischer
% Materialfließkurve definiert durch zyklische Ramberg-Osgood Parameter
%
% Funktion breitet Sützpunkte der RO Gleichung so auf, das diese der
% Funktion fk2para übergeben werden können. Sollte die Bauteilfließkurve
% berechnet werden so wird zusätzlich eine 1d Kerbnäherung durchgeführt um
% die Bauteilfließkurve zu bestimmen.
%
% _________________________________________________________________________
% INPUT:
%  typ            -> Welche Art parameter soll bestimmt werden (string)
%                       "werkstoff"
%                       "pseudo stress"
%                       "pseudo strain"
%                       "DevEps"
% E,nu            -> Konstanten isotrope Elastizität  (e R)
% Kstrich,nstrich -> Konstanten Ramberg Osgood gesetz (e R)
% M               -> Anzahl Backstresstensoren die verwendet werden sollen
%                                                     (e N)
% modell          -> verwendetes Materialmodell       (string)
%                      "Chaboche"
%                      "OhnoWang"
%                      "Jiang"
%                      "KarimOhno"
% verfahren_flag  -> Entscheidet, welches Verfahren zum setzten der
%                   Stützstellen verwendet wird.
% varargin        -> Variabler Input je nach Verfahren
%
% _________________________________________________________________________
%
% verschiedene Verfahren:
% verfahren_flag == 1 -> äquidistant in Spannunge
%                       varargin{1,1} = Rm -> Zugfestigkeit
% 
% verfahren_flag == 2 -> geometrische Folge nach Döring
%                       varargin{1,1} = q -> Quotien ep_i+1 = 1/q * ep_i
%
% verfahren_flag == 3 -> verfahren von Simon
%                       varargin{1,1} = q -> Anteil plastischer Dehnungen 
%                                            an gesamtdehnung bei
%                                            Fließbegindd
%                       varargin{1,2} = ep_M -> plastische dehunung nach
%                                               der ideale plastizität gilt
% _________________________________________________________________________
% 
% variabler Input (feste reihenfolge):
% q              -> Parameter zum Setzen der Stützstellen
% ep_M           -> Parameter zum Setzen der Stützstellen, bei
%                   verfahren_flag = 2 kann irgendein Wert gegeben werden
% ksim1d         -> Kerbnäherungsverfahren für 1d Kerbsimulation (string)
%                      "Neuber"
%                      "ESED"
%                      "Seeger Heuler"
%                      "Seeger Beste"
% Kp             -> Traglastformzahl für einige 1d Kerbnäherungen (e R)
%
% _________________________________________________________________________
% OUTPUT:
% para            -> Materialparameter des Modells
%
%
%
% -------------------------------------------------------------------------
% Stand: Januar 2021
% Autor: Jan Kraft
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Variablen Input 
% ... Erstmal default Werte
if verfahren_flag == 2 
    defaultval = {0.333, 0.03, 'Neuber', 1.001};
elseif verfahren_flag == 3
    defaultval = {0.01, 0.03, 'Neuber', 1.001};
elseif verfahren_flag == 1
    defaultval = {500, 0.03, 'Neuber', 1.001};
else
    verfahren_flag = 3;
    defaultval = {0.333, 0.03, 'Neuber', 1.001};
end
% ... setze variablen input
nvarin = size(varargin,2);
if nvarin > 0
    for i = 1:nvarin
        defaultval{i} = varargin{i};
    end
end
q = defaultval{1};
epM = defaultval{2};
ksim1d = defaultval{3}; 
Kp = defaultval{4}; 

% -------------------------------------------------------------------------
% initialisiere lokale variablen
sig = NaN(1,M+1);                                                          % spannungen
epsp = NaN(1,M+1);                                                         % plastische Dehnungen 


% -------------------------------------------------------------------------
% festlegen Ratchetting Parameter (default werte)
switch typ
    case 'werkstoff'
        chi_i = 5 * ones(1,M);                                             % Default werkstoff
    otherwise
        chi_i = 50 .* ones(1,M);                                           % Default struckturmodell
end

% -------------------------------------------------------------------------
% Äquidistante Stützstellen in Spannungen
if verfahren_flag == 1
    % Zug festigkeit
    Rm = q;%varargin{1,1};
    % zyklische Fließspannung definiert bei 0.01% plastischer Dehnung
    r0 = Kstrich * 0.0001^nstrich;

    % Spannungen und Dehnungen bei Fließbeginn
    sig(1) = r0;
    epsp(1) = 0.0001;

    % Stüzstellen echte zyklische Werkstoffkurve
    % angepasst an Zugfestigkeit
    for ii = 1 : M
        % stützstellen in spannungen (gleiche schrittweite)
        sig(ii+1) = r0 + (Rm - r0)/M*ii;
        % stützstellen plastische Dehnungen
        epsp(ii+1) = (sig(ii+1)/Kstrich)^(1/nstrich);
    end
    epsp(1) = 0.0;


% -------------------------------------------------------------------------
% Geometrische Folge
elseif verfahren_flag == 2
    % zyklische Fließspannung definiert bei 0.01% plastischer Dehnung
    r0 = Kstrich * 0.0001^nstrich;
    % Spannungen und Dehnungen bei Fließbeginn
    sig(1) = r0;
    epsp(1) = 0.0001;
    % Stüzstellen echte zyklische Werkstoffkurve
    % als geometrische Reihe nach Döring
    for ii = 1 : M
        % stützstelle plastische Dehnung
        epsp(ii+1) = epsp(ii)/q;
        % sützstelle spannungen
        sig(ii+1) = Kstrich * (epsp(ii+1))^nstrich;
    end

    epsp(1) = 0.0;
    
    
% -------------------------------------------------------------------------
% Verfahren von Simon
elseif verfahren_flag == 3
    % bestimme Fließspannung
    r0 = 10^( ( log10(q/((1-q)*E)) + 1/nstrich * log10(Kstrich) ) /... 
              (1/nstrich - 1) ...
              );
    epsp(1) = 0.0;
    sig(1) = r0;
    % bestimme ersten stützpunkt
    sig(2) = 1/(1-nstrich) * r0;
    epsp(2) = (sig(2)/Kstrich)^(1/nstrich);
    % bestimme Restliche Stützstellen
    for ii = 3 : M+1
        epsp(ii) = 10 ^ (log10(epsp(2)) + (ii-2) * ( log10(epM) - ...
                                                   log10(epsp(2)) )/(M-1));
        sig(ii) = Kstrich * epsp(ii)^nstrich;
    end
end


% -------------------------------------------------------------------------
% Pseudo Sützstellen
if ~strcmp(typ,'werkstoff')
    switch ksim1d
        case 'Neuber'
            [esig,eepsp]= uniax_neuber(E,sig,epsp);
        case 'ESED'
            [esig,eepsp]= uniax_esed(E,nstrich,sig,epsp);
        case 'Seeger Beste'
            [esig,eepsp] = uniax_seegerbeste(E,sig,epsp,Kstrich,nstrich,Kp);
        case 'Seeger Heuler'
            [esig,eepsp] = uniax_neuberstern(E,sig,epsp,Kstrich,nstrich,Kp);
        case 'Neuber Stern'
            [esig,eepsp] = uniax_neuberstern(E,sig,epsp,Kstrich,nstrich,Kp);
        otherwise
            msg = 'angegebenes Verfahren nicht implementiert';
            error(msg)
    end
end

% -------------------------------------------------------------------------
% bestimme Parameter aus Fließkurve
switch typ
    case 'werkstoff'
        para = fk2para(sig, epsp, E, nu, r0, chi_i, M, modell);
    case 'pseudo stress'
        para = fk2para(esig, epsp, E, nu, r0, chi_i, M, modell);
    case 'Pseudo Stress'
        para = fk2para(esig, epsp, E, nu, r0, chi_i, M, modell);
    case 'PseudoStress'
        para = fk2para(esig, epsp, E, nu, r0, chi_i, M, modell);
    case 'PseudoStress Lang'
        para = fk2para(esig, epsp, E, nu, r0, chi_i, M, modell);
    case 'pseudo strain'
        sig(end) = sig(end) * 0.99; % Abfangen numerischer fehler bei Spannungsintegration im pseudo strain approach
        para = fk2para(sig, eepsp, E, nu, r0, chi_i, M, modell);
    case 'Pseudo Strain'
        sig(end) = sig(end) * 0.99; % Abfangen numerischer fehler bei Spannungsintegration im pseudo strain approach
        para = fk2para(sig, eepsp, E, nu, r0, chi_i, M, modell);
    case 'PseudoStrain'
        sig(end) = sig(end) * 0.99; % Abfangen numerischer fehler bei Spannungsintegration im pseudo strain approach
        para = fk2para(sig, eepsp, E, nu, r0, chi_i, M, modell);
    case 'DevEps'
        G = E/(2*(1+nu));
        eepspv = epsp + (sig - esig)./(3*G);
        para = fk2para(esig, eepspv, E, nu, r0, chi_i, M, modell);
    otherwise
        msg = 'angegebenes Verfahren nicht implementiert';
        error(msg)
end

end % Ende Funktion




