function [ZVAR0] = init_matmodel(...
              matmodell, ntens, parameter, M)
% Funktion Initialisiert Materialmodel
% - Überprüft ob korrekte Anzahl der Materialparameter
% - Initialisiert Zustandsvariablen
% - Wählt Materialfunktion aus
%
% INPUT:
% matmodell    - (str) Name des Materialmodells
% ntens        - (int) Anzahl Tensorkomponenten
% parameter    - (array) Parameter des Materialmodells
% M            - (int) Anzahl der Backstresstensoren
%
% OUTPUT:
% ZVAR0        - Zustandsvariablen für unbelasteten Zustand
% -------------------------------------------------------------------------

switch matmodell
    case 'Chaboche'
        % Prüfe Parameter
        Msoll = (length(parameter)-5)/2;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter für Chaboche',...
                   ' übergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+2,1);
        r0 = parameter(end); %Anfangsradius der ff
        ZVAR0(end-1) = r0;
        
    case 'OhnoWang'
        % Prüfe Parameter
        Msoll = (length(parameter)-3)/3;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter für OhnoWang',...
                   ' übergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+1,1);    
        
    case 'KarimOhno'
        % Prüfe Parameter
        Msoll = (length(parameter)-3)/3;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter für KarimOhno',...
                   ' übergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+1,1); 
        
    case 'Jiang'
        % Prüfe Parameter
        Msoll = (length(parameter)-8)/7;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter für Jiang',...
                   ' übergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+2,1);
        
    case 'Doring'
        % Prüfe Parameter
        Msoll = (length(parameter)-22)/8;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter für Doring',...
                   ' übergeben'];
            error(msg)
        end
        
        % init Startwerte
        ZVAR0 = initZVARDoring(ntens,parameter,M);
    
    case 'OWT'
        % Prüfe Parameter
        Msoll = (length(parameter)-11)/3;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter für OWT',...
                   ' übergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((3+M+(ntens+1)/2)*ntens+4,1); 
        
    otherwise
        msg = 'Angegebenes Modell nicht implementiert';
        error(msg)
end
end