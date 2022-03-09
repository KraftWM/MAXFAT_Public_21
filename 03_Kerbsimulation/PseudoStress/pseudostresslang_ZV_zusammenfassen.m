function [ZVL] = pseudostresslang_ZV_zusammenfassen(ZVAR,EZVAR,material,M,eM)
% Funktion fasst Zustandsvariablen des Material- und Strukturmodells
% zusammen um den Pseudo Spannungsansatz mit dem Algo. nach Lang
% auszuführen
%
% INPUT:
%  ZVAR       - Zustandsvariablen des Materialmodells
%  EZVAR      - Zustandsvariablen des Strukturmodells
%  material   - Name des Plastizitätsmodells
%   M         - Anzahl der Backstresstensoren
%
% OUTPUT:
%  ZVL        - Zusammengefasste Zustandsvariablen
%
% _________________________________________________________________________

% !!!!!!!! nur ESZ
ntens = 3;


% =========================================================================
% Definiere Materialfunktion und Zustandsvariablen
% Zuerst Zustandsvariablen des Struckturmodells (wie bei Spannungssteuerung)
% dann zusätzliche Variablen des Materialmodells 
switch material
    
    case 'Chaboche'
        
        % Radius keine ZV -> Aus Lösung DGL
%         % Startwerte der Zustandsvariablen
%         ZVL(1:6+3*M,1) = EZVAR(1:6+3*M,1); % eeps,eepsp,ealphai
%         ZVL(7+3*M) = EZVAR(8+3*M);         % p
%         ZVL(8+3*M:7+6*M,1) = ZVAR(1+2*ntens:(2+M)*ntens,1); % alphai
        
        % Radius als ZV
        % Startwerte der Zustandsvariablen
        ZVL(1:6+3*eM,1) = EZVAR(1:6+3*eM,1); % eeps,eepsp,ealphai
        ZVL(7+3*eM,1) = EZVAR(7+3*eM,1); % eY
        ZVL(8+3*eM) = EZVAR(8+3*eM);         % p
        ZVL(9+3*eM:8+3*eM+3*M,1) = ZVAR(1+2*ntens:(2+M)*ntens,1); % alphai
        ZVL(9+3*eM+3*M,1) = ZVAR((2+M)*ntens+1,1);
        
    case 'OhnoWang'
        
        % Startwerte der Zustandsvariablen
        ZVL(1:7+3*eM,1) = EZVAR;
        ZVL(8+3*eM:7+3*eM+3*M,1) = ZVAR(1+2*ntens:(2+M)*ntens,1);

        
    case 'KarimOhno'
        
        % Startwerte der Zustandsvariablen
        ZVL(1:7+3*eM,1) = EZVAR;
        ZVL(8+3*eM:7+3*eM+3*M,1) = ZVAR(1+2*ntens:(2+M)*ntens,1);

        
    case 'Jiang'
        
        % Startwerte der Zustandsvariablen
        ZVL(1:8+3*eM,1) = EZVAR;
        ZVL(9+3*eM:8+3*eM+3*M,1) = ZVAR(1+2*ntens:(2+M)*ntens,1);
        ZVL(9+3*eM+3*M,1) = ZVAR(end,1);

    case 'OWT'
        
        % Startwerte der Zustandsvariablen
        ZVL(1:(5+eM)*3+4,1) = EZVAR;
        ZVL((5+eM)*3+5:(8+M+eM)*3+7,1) = ZVAR([2*ntens+1:(2+M)*ntens , ...                    % alphai
                                             (2+M)*ntens+2 , ...                            % Q
                                              (2+M)*ntens+3:(3+M)*ntens+2, ...              % beta
                                              (3+M)*ntens+3, ...                            % q
                                              (3+M)*ntens+4, ...                            % A
                                              (3+M)*ntens+5 : (3+M+(ntens+1)/2)*ntens+4 ... % CT
                                              ],1);
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert für Ansatz nach Lang';
        error(msg)
end
