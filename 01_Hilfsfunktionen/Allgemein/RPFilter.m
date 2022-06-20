function [L,DLZ,verk] = RPFilter(L,level,nkomb)
% Rainflow Projection Filter nach:
%
% QUELLE:
%
% 
% INPUT:
% L      - Lastsignal (nkana x ndata)
% level  - relatives Filter Level
% nkomb  - Anzahl an Linearkombinationen
%
% OUTPUT:
% L      - Verkuerztes Lastsignal
% verk   - Verkuerzung des Signals 
%
% -------------------------------------------------------------------------

%% Eigenschaften des Signal
nkana = size(L,1);
ndata = size(L,2);
RES = zeros(ndata,1);
DLZ = 1:ndata;

%% Erzeuge Linearkombinationen
% c = linearKombinationen(nkana,nkomb);
if nkana == 2 && nkomb ~= 1
    nkomb = 9;
    c =  [ 1.00, 0.00;...
     	   0.89, 0.45;... 
     	   0.71, 0.71;...
     	   0.45, 0.89;... 
     	   0.00, 1.00;... 
     	  -0.45, 0.89;... 
     	  -0.711, 0.709;... % geringe Verstimmung fuer proportionale Belastung notwendig
     	  -0.89, 0.45;...
     	  -1.00, 0.00];
else 
    nkomb = 1;
    c = 1;
end

% Zu behaltende Werte
keep = nkomb * ones(ndata,1);

%% Kombinierte Signale
% Schleife ueber kombinierte Signale
for i = 1 : nkomb
    Lkomb = c(i,:) * L;                                                    % Kombiniertes Lastsignal
    Lamax = (max(Lkomb)-min(Lkomb))/2;                                     % Maximale Amplitude
%     AMP = hcm_origV2(Lkomb);
%     Lamax = max(AMP(5,:));
    RES(1) = 1; IZ = 1; IR = 1; counter = 0;                               % Init Rainflow
    % Rainflow
    for K = 1:ndata
        [IZ,IR,RES,counter,keep] = hcm(...
                              K,IZ,IR,RES,counter,Lkomb,keep,Lamax,level);
    end
    % Hoechstwert nochmal aufbringen
    [~,~,RES,~,keep] = hcm(...
                              RES(IR),IZ,IR,RES,counter,Lkomb,keep,Lamax,level);
end

%% Signal Verkuerzen
keep = keep >= 1;
L = L(:,keep);
DLZ = DLZ(keep);
verk = size(L,2)/ndata;



end % Ende Hauptfunktion


% -------------------------------------------------------------------------
% Hilfsfunktion erzeuge Linearkombinationen
function c = linearKombinationen(n,k)
% !!!!! Aktuell nur fuer n = 1,2
if n > 2
    error('RPFilter aktuell nur bis n=2')
end
% Init Koeffizienten
N = (k+1)/2 * k^(n-1) - (k-1)/2 * (k-2)^(n-1);
c = zeros(N,n);
% Uniformes Netz

% Normiere Netz


end % Ende Erzeuge Linearkombinationen



% -------------------------------------------------------------------------
%   Hilfsfunktion HCM fuer Filterung                 
function [IZ,IR,RES,counter,keep] = hcm(...
                               K,IZ,IR,RES,counter,Data,keep,Lamax,level)

    % Toleranzen
    tolM12 = 0.99;                             % Toleranz fuer das Erkennen vom Memory 1 und 2
    tolM3 = 1.01;                              % Toleranz fuer das Erkennen vom Memory 3

    % Abbruchbedingung
    weiter = 1;

    % ---------------------------------------------------------------------
    % Rainflow
    while weiter
        % 2
        if IZ > IR % Vergleich der Zeiger 
            
            % ... letzte Werte aus Residuum lesen 
            I = RES(IZ-1);
            J = RES(IZ);
            
            % ... Pruefe ob letzter Wert UKP ist
            if (Data(K)-Data(J))*(Data(J)-Data(I)) >= 0
                % ... kein UKP
                IZ = IZ - 1;
                weiter = 1;
%                 !! GOTO 2
            else
                % ... UKP
                % ... Pruefe Schwingweite groeßer als die letzte
                if abs(Data(K)-Data(J)) >= tolM12 * abs(Data(J)-Data(I))
                    % ... Schwingspiel gefunden
                    counter = counter + 1;
                    AMP = 0.5* abs(Data(I)-Data(J));
                    % Amplitude kleiner als Filterlevel ?
                    if AMP < Lamax * level
                        % Hystrese Ausschneiden
%                         keep(I+1:K-1) = keep(I+1:K-1) - 1;
                        keep(I:K-1) = keep(I:K-1) - 1;
                    end
                    IZ = IZ - 2;
                    % ... Stapel leer ?
                    if IZ >= IR
                        % ... nein
%                         !! GOTO 2
                        weiter = 1;
                    else
                        % ... ja
                        IZ = IZ + 1;
                        weiter = 0;
                    end
                else
                    % ... kein Schwingspiel
                    weiter = 0;
                    IZ = IZ + 1;
                end % Ende Verzweigung ueberpruefung der Schwingweiten
            end % Ende Verzweigung UKP
            
        else
            % ... IZ <= IR (wird hier zsm behandelt)
            % ... Einlesen Wert aus Residuum
            J = RES(IZ);
            weiter = 0;
            % ... Pruefe UKP
            if (Data(K)-Data(J))*Data(J) < 0
                % ... UKP
                % ... Pruefe Memory 3
                if abs(Data(K)) > tolM3 * abs(Data(J))
                    % ... Memory 3
                    IR = IR + 1;
                end
                IZ = IZ + 1;
            end
        end % Ende Verzweigung Zeigervergleich
    end % Ende while Schleife
    
    RES(IZ) = K;
    
end % Ende HCM Funktion











