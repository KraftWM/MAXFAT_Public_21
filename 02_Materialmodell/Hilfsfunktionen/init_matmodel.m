function [ZVAR0] = init_matmodel(...
              matmodell, ntens, parameter, M)
% Funktion Initialisiert Materialmodel
% - �berpr�ft ob korrekte Anzahl der Materialparameter
% - Initialisiert Zustandsvariablen
% - W�hlt Materialfunktion aus
%
% INPUT:
% matmodell    - (str) Name des Materialmodells
% ntens        - (int) Anzahl Tensorkomponenten
% parameter    - (array) Parameter des Materialmodells
% M            - (int) Anzahl der Backstresstensoren
%
% OUTPUT:
% ZVAR0        - Zustandsvariablen f�r unbelasteten Zustand
% -------------------------------------------------------------------------

switch matmodell
    case 'Chaboche'
        % Pr�fe Parameter
        Msoll = (length(parameter)-5)/2;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter f�r Chaboche',...
                   ' �bergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+2,1);
        r0 = parameter(end); %Anfangsradius der ff
        ZVAR0(end-1) = r0;
        
    case 'OhnoWang'
        % Pr�fe Parameter
        Msoll = (length(parameter)-3)/3;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter f�r OhnoWang',...
                   ' �bergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+1,1);    
        
    case 'KarimOhno'
        % Pr�fe Parameter
        Msoll = (length(parameter)-3)/3;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter f�r KarimOhno',...
                   ' �bergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+1,1); 
        
    case 'Jiang'
        % Pr�fe Parameter
        Msoll = (length(parameter)-8)/7;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter f�r Jiang',...
                   ' �bergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((2+M)*ntens+2,1);
        
    case 'Doring'
        % Pr�fe Parameter
        Msoll = (length(parameter)-22)/8;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter f�r Doring',...
                   ' �bergeben'];
            error(msg)
        end
        
        % init Startwerte
        ZVAR0 = initZVARDoring(ntens,parameter,M);
    
    case 'OWT'
        % Pr�fe Parameter
        Msoll = (length(parameter)-11)/3;
        if M ~= Msoll
            msg = ['Falsche Anzahl an Materrialparameter f�r OWT',...
                   ' �bergeben'];
            error(msg)
        end
        % init Startwerte
        ZVAR0 = zeros((3+M+(ntens+1)/2)*ntens+4,1); 
        
    otherwise
        msg = 'Angegebenes Modell nicht implementiert';
        error(msg)
end
end