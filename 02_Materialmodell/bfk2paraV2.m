function para = bfk2paraV2(typ,bfk,M,modell,E,nu,verfahren_flag,varargin)
% Funktion zum bestimmen der Materialparameter/ Struckturparameter aus
% Bauteilfließkurve.
%
% funktion bereitet bfk so auf, dass diese der funktion fk2para übergeben
% werden kann
%
% Dafür das die Fließkurve zum jeweiligen ansatz passt muss der Anwender
% selber sorgen
%
% _________________________________________________________________________
% INPUT:
%  typ            -> Welche Art parameter soll bestimmt werden (string)
%                       "werkstoff"
%                       "pseudo stress"
%                       "pseudo strain"
%                       "DevEps"
% M               -> Anzahl Backstresstensoren die verwendet werden sollen
%                                                     (e N)
% bfk             -> Bauteilfließkurve (e R^(numdat x 2))
%                 1. Spalte = Spannungen   2. Spalte = plastische Dehnungen   
% modell          -> verwendetes Materialmodell       (string)
%                      "Chaboche"
%                      "OhnoWang"
%                      "Jiang"
%                      "KarimOhno"
% E, nu           -> Elastizitätsparameter
% verfahren_flag -> Entscheidet, welches Verfahren zum setzten der
%                   Stützstellen verwendet wird.
%
% _________________________________________________________________________
% OUTPUT:
% para            -> Materialparameter des Modells
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
%                       varargin{1,3} = Kstrich -> RamOsg.
%                       varargin{1,4} = nstrich -> RamOsg.
% -------------------------------------------------------------------------
% Stand: September 2020
% Autor: Jan Kraft
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% herrauslesen Spannungen und plastische Dehnungen
sig = bfk(:,1);
epsp = bfk(:,2);
numdat = size(sig,1);

% -------------------------------------------------------------------------
% festlegen Ratchetting Parameter (default werte)
switch typ
    case 'werkstoff'
        chi_i = 5 * ones(1,M);                                             % Default werkstoff
    otherwise
        chi_i = 50 .* ones(1,M);                                           % Default struckturmodell
end


% -------------------------------------------------------------------------
% Stützstellen für Parameter Ermittlung
if verfahren_flag == 1 
    % .. äquidistante Stützstellen 
    
    % zyklische Fließspannung definiert bei 0.01% plastischer Dehnung
    % Festlegend Fließspannung bei ersten wert an dem epsp < 0.0001
    epspf = 0.0001; % plastische Dehnung an der Fließen definiert wird
    weiter = 1;     % Bool bedingung
    iter = 1;       % schleifen zähler
    while weiter
        % keine Endlosschleifen
        if iter > numdat - 1
            msg = ['zu kleine plastischen Dehnngen in Bauteilfließkurve!',...
                   ' Parameter können nicht bestimmt werden'];
            error(msg)
        end
        % festlegen Fließspannung wenn epspf überschritten wird
        if epsp(iter) <= epspf && epsp(iter+1) > epspf
            % Abbruchbedingung
            weiter = 0;
            % Fließspannung
            r0 = sig(iter);
            % Zurücksetzten für späteren zugriff
            iter = iter - 1;
        elseif iter == 1 && epsp(iter) > epspf
            % Fließspannung
            r0 = sig(1);
            % Warnung ausgeben
            msg = ['Erste plastische Dehnung in Bauteilfließkurve größer als ',...
                   'festgelegte Grenze (epsp_f = ',num2str(epspf),') für' ,...
                   'Fließen. Erster Spannungswert wird ', ...
                   'als Fließspannung übernommen. epsp(1) = ', num2str(epsp(1)), ...
                   ' Fließspannung sig_f = ', num2str(r0)];
            warning(msg);
            % Abbruchbedingung
            weiter = 0;
        end
        iter = iter + 1;
    end
    
    % Grenzwert
    Rm = varargin{1,1};

    % Spannungen und Dehnungen bei Fließbeginn
    sigwerte = zeros(1,M+1);                 % Stützstellen Spannungen
    sigwerte(1) = r0;
    epspwerte = zeros(1,M+1);                % Stützstellen plastische Dehnungen
    epspwerte(1) = epsp(iter);

    % Stüzstellen echte zyklische Werkstoffkurve
    % angepasst an Zugfestigkeit
    for ii = 1 : M
        % stützstellen in spannungen (gleiche schrittweite)
        signext = r0 + (Rm - r0)/M*ii;
        % finde
        idx = find(sig > signext,1,'first');    % erster Wert in sig >= r0
        if isempty(idx) % Gesuchter Wert Nicht in Lastfolge
            % ... Ausgabe Warnung
            msg = ['Stützpunkt ', num2str(ii+1), ' wird extrapoliert. ',...
                   'Aktueller Spannungswert ', num2str(signext), ' höher als ',...
                   'maximaler Wert in gegebener Fließkurve ',num2str(max(sig))];
            warning(msg);
            
            % Stützstelle Spannungen
            sigwerte(ii+1) = signext;
            % stützstellen plastische Dehnungen mit extrapolation aus
            % letzten beiden Werten
            epspwerte(ii+1) = linint(sig,epsp,length(sig),signext);
        else
            % Stützstelle Spannungen
            sigwerte(ii+1) = signext;%sig(idx);
            % stützstellen plastische Dehnungen
            epspwerte(ii+1) = linint(sig,epsp,idx,signext);%epsp(idx);    % erster Wert in sig >= r0;
        end
        
    end
    epspwerte(1) = 0.0;
    
elseif verfahren_flag == 2
    % ... stützstellen als geometrische Folge
    
    % Festlegend Fließspannung bei ersten wert an dem epsp < 0.0001
    epspf = 0.0001; % plastische Dehnung an der Fließen definiert wird
    weiter = 1;     % Bool bedingung
    iter = 1;       % schleifen zähler
    while weiter
        % keine Endlosschleifen
        if iter > numdat - 1
            msg = ['zu kleine plastischen Dehnngen in Bauteilfließkurve!',...
                   ' Parameter können nicht bestimmt werden'];
            error(msg)
        end
        % festlegen Fließspannung wenn epspf überschritten wird
        if epsp(iter) <= epspf && epsp(iter+1) > epspf
            % Abbruchbedingung
            weiter = 0;
            % Fließspannung
            r0 = sig(iter);
            % Zurücksetzten für späteren zugriff
            iter = iter - 1;
        elseif iter == 1 && epsp(iter) > epspf
            % Fließspannung
            r0 = sig(1);
            % Warnung ausgeben
            msg = ['Erste plastische Dehnung in Bauteilfließkurve größer als ',...
                   'festgelegte Grenze (epsp_f = ',num2str(epspf),') für' ,...
                   'Fließen. Erster Spannungswert wird ', ...
                   'als Fließspannung übernommen. epsp(1) = ', num2str(epsp(1)), ...
                   ' Fließspannung sig_f = ', num2str(r0)];
            warning(msg);
            % Abbruchbedingung
            weiter = 0;
        end
        iter = iter + 1;
    end

    sigwerte = zeros(1,M+1);                 % Stützstellen Spannungen
    sigwerte(1) = r0;
    epspwerte = zeros(1,M+1);                % Stützstellen plastische Dehnungen
    epspwerte(1) = epsp(iter);
    q = varargin{1,1};                       % Quotient für geometrische Reihe
    epsp0 = epsp(iter);                      % Alte Stützstelle
    % Schleife für alle Backstresstensoren
    for ii = 1 : M
        % Neue Stützstelle nach geometrischer Reihe
        epsp1 = epsp0/q;
        % Finde plastischen Dehnungswert, der am nächsten an epsp1 ist
        idx = find(epsp > epsp1,1,'first');    % erster Wert in epsp >= epsp1
        if isempty(idx) % Gesuchter Wert Nicht in Lastfolge
            % ... Ausgabe Warnung
            msg = ['Stützpunkt ', num2str(ii+1), ' wird extrapoliert. ',...
                   'Aktueller Dehnungswert ', num2str(epsp1), ' höher als ',...
                   'maximaler Wert in gegebener Fließkurve ',num2str(max(epsp))];
            warning(msg);
            
            % Stützstelle Spannungen
            epspwerte(ii+1) = epsp1;
            % stützstellen plastische Dehnungen mit extrapolation aus
            sigwerte(ii+1) = linint(epsp,sig,length(sig),epsp1);
        else
            % stützstellen plastische Dehnungen
            epspwerte(ii+1) = epsp1;%epsp(idx); 
            % Stützstelle Spannungen
            sigwerte(ii+1) = linint(epsp,sig,idx,epsp1);%sig(idx);              
        end
        % speichern plastische Stützstelle
        epsp0 = epsp1;

    end 
    epspwerte(1) = 0;
elseif verfahren_flag == 3
    % ... Verfahren von Simon
    % plastischer Dehnunsanteil
    q = varargin{1,1};
    % letzter Stützpunkt plastische Dehnungen
    epM = varargin{1,2};
    % Ramberg Osg. Parameter 
    idx = find(sig < (1-q)*(sig + epsp*E),1,'first'); 
    if isempty(idx)
        msg = ['zu große Spannungen in Bauteilfließkurve!',...
                   ' Parameter können nicht bestimmt werden'];
            error(msg)
    else
        % init andere Stützstellen
        r0 = sig(idx);
        sigwerte = zeros(1,M+1);                 % Stützstellen Spannungen
        sigwerte(1) = sig(idx);
        epspwerte = zeros(1,M+1);                % Stützstellen plastische Dehnungen
        epspwerte(1) = epsp(idx);
    end
    
    % finde andere Stützstellen (Punkt dessen Tangente durch die
    % Fließspannungen geht)
    H = (sig(2:end)-sig(1:end-1))./(epsp(2:end)-epsp(1:end-1));
    idx2 = find( r0 < -H.*epsp(1:end-1)+sig(1:end-1),1,'first'); 

    if idx2 + 1 < idx || isempty(idx2)
        msg = 'etwas ist schief gelaufen';
        error(msg)
    else
        sigwerte(2) = sig(idx2);
        epspwerte(2) = epsp(idx2);
    end
    
    
    % Finde restliche Stützpunkte
    for ii = 3:M+1
        epsp1 = 10 ^ (log10(epspwerte(2)) + (ii-2) * ( log10(epM) - ...
                                              log10(epspwerte(2)) )/(M-1));
        idx = find(epsp > epsp1,1,'first');
        if isempty(idx)
            msg = ['Stützpunkt ', num2str(ii+1), ' wird extrapoliert. ',...
                   'Aktueller Dehnungswert ', num2str(epsp1), ' höher als ',...
                   'maximaler Wert in gegebener Fließkurve ',num2str(max(epsp))];
            warning(msg)
            
            epspwerte(ii) = epsp1;
            sigwerte(ii) = linint(epsp,sig,length(epsp),epsp1);
        else
            epspwerte(ii) = epsp1;
            sigwerte(ii) = linint(epsp,sig,idx,epsp1);
        end
            
    end
    
    epspwerte(1) = 0;
end

% -------------------------------------------------------------------------
% bestimme Parameter aus Fließkurve
para = fk2para(sigwerte, epspwerte, E, nu, r0, chi_i, M, modell);


end % Ende Hauptfunktion 


% Hilfsfunktion lineare interpolation
function yw = linint(x,y,idx,xw)
x1 = x(idx-1);
x2 = x(idx);
y1 = y(idx-1);
y2 = y(idx);
yw = (y2-y1)/(x2-x1)*(xw-x1)+y1;
end
