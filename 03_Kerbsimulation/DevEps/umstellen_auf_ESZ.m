function [ZVARESZ] = umstellen_auf_ESZ(ZVAR,material,M)
% Funktion erhält Zustandsvariablen im 3D und stellt sie um auf ESZ zustand
% (nötig für DevEps ansatz)
%
% (dehnungsgesteuertes) Werkstoffmodell wird auf ESZ Modell
% umgestellt
% (spannungsgesteuertes) Strukturmodell wird auf ESZ Modell
% umgestellt
%
% INPUT:
% ZVAR             - Zustandsvariablen Materialmodell im 3D
% material         - Name des Plastizitätsmodells
% M                - Anzahl der Backstresstensoren 
% flag             - 0 Spannungssteuerung 1 Dehnungssteurung
%
% OUTPUT:
% ZVARESZ          - Zustandsvariablen für ESZ
%
%__________________________________________________________________________


% neue Tensoranzahl Speichern
ntens = 3;
% anzahl Datenpunkte in Zeitreihe
numink = size(ZVAR,2);

% Anzahl Zustandsvariablen
switch material
    case 'Chaboche'
        numzvar = (2+M)*3+2;
    case 'KarimOhno'
        numzvar = (2+M)*3+1;
    case 'OhnoWang'
        numzvar = (2+M)*3+1;
    case 'Jiang'
        numzvar = (2+M)*3+2;
    case 'OWT'
        numzvar = (3+M+(3+1)/2)*3+4;
    case 'Doring'
        msg = 'Döring Modell hier noch net implementiert';
        error(msg)       
    otherwise        
        msg = ['Angegebenes Materialmodell, ',material,', nicht implementiert'];
        error(msg)
end

% Speicher Zustandsvariablen
ZVARESZ = zeros(numzvar,numink);

% Spannungen/Dehnungen
ZVARESZ(1:ntens,:) = ZVAR([1,2,4],:);

% plastische Dehnungen
ZVARESZ(ntens+1:2*ntens,:) = ZVAR([7,8,10],:);

% Backstress
for jj = 2 : M + 1
    idx3D = jj*6+1;
    ZVARESZ(jj*ntens+1:(jj+1)*ntens,:) = ZVAR([idx3D, idx3D+1, idx3D+3],:);
end
        
% materialfunktion
switch material
    
    % Definiere Material Funktion
    case 'Chaboche'             
        % radius FF
        ZVARESZ(end-1,:) = ZVAR(end-1,:);
        
        % plastische Bogenlänge
        ZVARESZ(end,:) = ZVAR(end,:);
                
    case 'OhnoWang'        
        % plastische Bogenlänge
        ZVARESZ(end,:) = ZVAR(end,:);
    
    case 'KarimOhno'
        % plastische Bogenlänge
        ZVARESZ(end,:) = ZVAR(end,:);
        
    case 'Jiang'        
        % plastische Bogenlänge
        ZVARESZ(end-1,:) = ZVAR(end-1,:);
        
        % radius gedächtnissfläche
        ZVARESZ(end,:) = ZVAR(end,:);
        
    case 'Doring'
        
        msg = 'Döring Modell ist für den DevEps Ansatz noch nicht implementiert';
        error(msg)
        
    case 'OWT'
        
        % plastische Bogenlänge
        ZVARESZ((2+M)*ntens+1,:) = ZVAR((2+M)*6+1,:);
        
        % Zusätzlicher FF Radius
        ZVARESZ((2+M)*ntens+2,:) = ZVAR((2+M)*6+2,:);
        
        % Backstraintensor
        idx3D = (2+M)*6+3;
        ZVARESZ((2+M)*ntens+3:(3+M)*ntens+2,:) = ZVAR([idx3D, idx3D+1, idx3D+3],:);
        
        % Radius gedächtnissfläche
        ZVARESZ((3+M)*ntens+3,:) = ZVAR((3+M)*6+3,:);
        
        % NP Parameter
        ZVARESZ((3+M)*ntens+4,:) = ZVAR((3+M)*6+4,:);
        
        % Tanaka Tensor
        idx = (3+M)*6+5;
        idx3D = [idx idx+1 idx+3 idx+6 idx+12 idx+15];
        ZVARESZ((3+M)*ntens+5 : (3+M+(ntens+1)/2)*ntens+4,:) = ZVAR(idx3D,:);
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
end

end % ENde Funktion