function [L,DLZ,verk] = RPMinMaxFilter(L,nkomb)
% Rainflow Projection Umkehrpunkt Filter nach:
%
% QUELLE:
%
% 
% INPUT:
% L      - Lastsignal (nkana x ndata)
% nkomb  - Anzahl an Linearkombinationen
%
% OUTPUT:
% L      - Verkürztes Lastsignal
% verk   - Verkürzung des Signals 
%
% -------------------------------------------------------------------------

%% Eigenschaften des Signal
nkana = size(L,1);
ndata = size(L,2);
DLZ = 1:ndata;

%% Erzeuge Linearkombinationen
% c = linearKombinationen(nkana,nkomb);
if nkana == 2
    nkomb = 9;
    c =  [ 1.00, 0.00;...
     	   0.89, 0.45;... 
     	   0.71, 0.71;...
     	   0.45, 0.89;... 
     	   0.00, 1.00;... 
     	  -0.45, 0.89;... 
     	  -0.711, 0.709;... % geringe Verstimmung für proportionale Belastung notwendig
     	  -0.89, 0.45;...
     	  -1.00, 0.00];
else 
    nkomb = 1;
    c = 1;
end

% Zu behaltende Werte
keep = false(1,ndata);

%% Kombinierte Signale
% Schleife über kombinierte Signale
for i = 1 : nkomb
    Lkomb = c(i,:) * L;                                                    % Kombiniertes Lastsignal
    keep = myfilterUKP(Lkomb,keep);                                          % Filtern von UKP
end

%% Signal Verkürzen
L = L(:,keep);
DLZ = DLZ(keep);
verk = size(L,2)/ndata;
end % Ende Hauptfunktion


% -------------------------------------------------------------------------
% Hilfsfunktion filtern Umkehrpunkte
function keep = myfilterUKP(L,keep)
% Bilde 1. Ableitung
dL = diff(L,[],2);
% Vorzeichenwechsel
UKP = dL(:,1:end-1).*dL(:,2:end) < 0;
UKP = [true,UKP,true]; % Ersten und letzten Punkt immer behalten
% Speichern UKP
keep = keep | UKP;
end % Ende UKP Filter











