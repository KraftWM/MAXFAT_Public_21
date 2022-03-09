function [ZVAR3D] = umstellen_auf_3D(ZVAR,material,M,nu,flag)
% Funktion erhält Zustandsvariablen im ESZ und stellt sie um auf 3D zustand
% (nötig für DevEps ansatz)
%
% (dehnungsgesteuertes) Werkstoffmodell wird auf 3D Modell
% umgestellt (nur fürn ESZ)
% (spannungsgesteuertes) Struckturmodell wird auf 3D Modell
% umgestellt (nur fürn ESZ)
%
% INPUT:
% ZVAR             - Zustandsvariablen Materialmodell im ESZ
% material         - Name des Plastizitätsmodells
% M                - Anzahl der Backstresstensoren  
% nu               - Querdehnzahl
% flag             - 0 Spannungssteuerung 1 Dehnungssteurung
%
% OUTPUT:
% ZVAR3D            - Zustandsvariablen für 3D
%
%__________________________________________________________________________


% alte Tensoranzahl Speichern
ntens = 3;
% anzahl Datenpunkte in Zeitreihe
numink = size(ZVAR,2);

% Anzahl Zustandsvariablen
switch material
    case 'Chaboche'
        numzvar = (2+M)*6+2;
    case 'KarimOhno'
        numzvar = (2+M)*6+1;
    case 'OhnoWang'
        numzvar = (2+M)*6+1;
    case 'Jiang'
        numzvar = (2+M)*6+2;
    case 'OWT'
        numzvar = (3+M+(6+1)/2)*6+4;
    case 'Doring'
        msg = 'Döring Modell hier noch net implementiert';
        error(msg)       
    otherwise        
        msg = ['Angegebenes Materialmodell, ',material,', nicht implementiert'];
        error(msg)
end

% Speicher
ZVAR3D = zeros(numzvar,numink);

% Spannungen & Dehnungen
EPSP = ZVAR(ntens+1:2*ntens,:);                      % plast. Dehnungen 

% ... Unterscheide Spannungs & Dehnungssteuerung

if flag % Dehnugssteuerung
   ezzp = - EPSP(1,:) - EPSP(2,:);
   EPSP = [EPSP(1:2,:);ezzp;EPSP(3,:)];
   ZVAR3D([1,2,4],:) = ZVAR(1:ntens,:);
   ZVAR3D([7,8,9,10],:) = EPSP;
else % Spannungssteuerung
   EPS = ZVAR(1:ntens,:); 
   [EPS, EPSP] = dehnungZZ(EPS,EPSP,nu);
   ZVAR3D([1,2,3,4],:) =  EPS;              
   ZVAR3D([7,8,9,10],:) = EPSP;
end

% Backstress
for jj = 2 : M+1
    AI = ZVAR(jj*ntens+1:(jj+1)*ntens,:);
    aizz = -AI(1,:)-AI(2,:);
    AI = [AI(1:2,:);aizz;AI(3,:)];
    idx = jj * 6 + 1;
    ZVAR3D([idx idx+1 idx+2 idx+3],:) = AI;
end
        
% unterscheide materialfunktion
switch material
    
    % Definiere Material Funktion
    case 'Chaboche'      
        % radius FF
        ZVAR3D(end-1,:) = ZVAR(end-1);
        
        % plastische Bogenlänge
        ZVAR3D(end,:) = ZVAR(end,:);               
        
    case 'OhnoWang'        
        % plastische Bogenlänge
        ZVAR3D(end,:) = ZVAR(end,:);
    
    case 'KarimOhno'
        % plastische Bogenlänge
        ZVAR3D(end,:) = ZVAR(end,:);
        
    case 'Jiang'       
        % plastische Bogenlänge
        ZVAR3D(end-1,:) = ZVAR(end-1);
        
        % radius gedächtnissfläche
        ZVAR3D(end,:) = ZVAR(end,:);
        
    case 'OWT'
                
        % plastische Bogenlänge
        ZVAR3D((2+M)*6+1,:) = ZVAR((2+M)*ntens+1,:);
        
        % Zusätzlicher FF Radius
        ZVAR3D((2+M)*6+2,:) = ZVAR((2+M)*ntens+2,:);
        
        % Backstraintensor
        idx = (2+M)*6 + 3;
        B = ZVAR((2+M)*ntens+3:(3+M)*ntens+2,:);
        bzz = -B(1,:) - B(2,:);
        ZVAR3D([idx idx+1 idx+2 idx+3],:) = [B(1,:);B(2,:);bzz;B(3,:)];
        
        % Radius gedächtnissfläche
        ZVAR3D((3+M)*6+3,:) = ZVAR((3+M)*ntens+3,:);
        
        % NP Parameter
        ZVAR3D((3+M)*6+4,:) = ZVAR((3+M)*ntens+4,:);
        
        % Tanaka Tensor
        idx = (3+M)*6+5;
        idx3D = [idx idx+1 idx+3 idx+6 idx+12 idx+15];
        ZVAR3D(idx3D,:) = ZVAR((3+M)*ntens+5 : (3+M+(ntens+1)/2)*ntens+4,:);
                
    case 'Doring'
        
        msg = 'Döring Modell hier noch net implementiert';
        error(msg)        
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
end

end % ENde Funktion