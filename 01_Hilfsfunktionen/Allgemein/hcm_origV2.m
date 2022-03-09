function [AMP,HYMAT] = hcm_origV2(Data)
% HCM Hysteresis Counting Method um Amplituden in einem Signal zu finden
% 
% Im Residuum werden keine Werte sondern nur Spaltenindices
%       gespeichert 
%
% Der vorliegende Code wurde nach
% RAINFLOW-HCM Ein Hyseresisschleifen-Z�hlalgorithmus auf
% Werkstoffmechanischer Grundlage von U.H. Chlormann und T. Seeger aus dem
% Jahr 1985 implementiert
%
%=========================================================================
% INPUT:
%
%
% Aus Kerbsimulation
% Data              - Diskrete Lastpunkt f�r Spannungen/Dehungen/Lasten 
%                     
%
%
%
% OUTPUT:
% AMP                - Amplituden
%                      (1.Zeile) -> Schwingspielz�hler
%                      (2.Zeile) -> Index in Data an dem Hyst schlie�t
%                      (3.Zeile) -> Amplitude K-J
%                      (4.Zeile) -> Mittelwert K-J
%                      (5.Zeile) -> Amplitude J-I
%                      (6.Zeile) -> Mittelwert J-I
% HYMAT              - RAINFLOW Matrix (!!!! noch machen !!!!)
%=========================================================================
% Erstellt von Jan Kraft
% Version 2.0 September 2020


% -------------------------------------------------------------------------
% Bereitstellen Speicher 
HYMAT = 0;
Ndata = size(Data,2);                                                      % Anzahl Datenpunkte
RES = zeros(1,Ndata);                                                      % reservierter Speicher f�r Residuum Mode I
AMP = zeros(6,Ndata);                                                      % Sch�digungsparameter (hier Amplituden da das nur ne testfunktion is)


% -------------------------------------------------------------------------
% Initialisierungen
% ... Init Residuum
RES(1) = 1;                                                                % Startwert Residuum direkt auf 1
% ... Init Z�hler zum Speicher von Amplituden
pz = 0;                                                                    % zeiger auf n�chste frei Spalte in Pz
% ... Init Nachfolger
NF = 2:Ndata+1;
NF(Ndata) = 1;
% ... Init Vorg�nger
VG = 0:Ndata-1;
VG(1) = Ndata;
% ... Init Zeiger auf n�chste einzulesende Stelle in Data
nextData = 1;
% ... Zyklusz�hler
counter = 0;                                 % Zyklusz�hler
% ... Zeiger auf Residuum
IR = 1;                                    % Zeiger auf Residuum
IZ = 1;                                    % Zeiger auf Residuum


% -------------------------------------------------------------------------
% Rainflow Z�hlunf
% ... Startindex
K = NF(RES(IZ));
% ... Schleife �ber alle Werte
while K ~= nextData
    
    % ... HCM aufrufen
    [IZ,IR,RES,counter,AMP,pz,NF,VG] = hcm(K,IZ,IR,RES,counter,Data,AMP,pz,NF,VG);
    
    % ... Inkrementiere K
    K = NF(K);
end


AMP = AMP(:,AMP(1,:) ~= 0);  % Aussortieren leere PZ Werte

end % Ende Hauptfunktion
























% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
%                                                                         %
%                   Hilfsfunktion HCM f�r einzelne Werte                  %
%     -> Berechnet PZ f�r die einzelnen Moden                             %
%     -> aktualisiert Rissschlie�dehnungen etc                            %
%                                                                         %
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %

function [IZ,IR,RES,counter,AMP,pz,NF,VG] = hcm(K,IZ,IR,RES,counter,Data,AMP,pz,NF,VG)

    % Toleranzen
    tolM12 = 0.99;                             % Toleranz f�r das Erkennen vom Memory 1 und 2
    tolM3 = 1.01;                              % Toleranz f�r das Erkennen vom Memory 3

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
            
            % ... Pr�fe ob letzter Wert UKP ist
            if (Data(K)-Data(J))*(Data(J)-Data(I)) >= 0
                % ... kein UKP
                IZ = IZ - 1;
                weiter = 1;
%                 !! GOTO 2
            else
                % ... UKP
                % ... Pr�fe Schwingweite gr��er als die letzte
                if abs(Data(K)-Data(J)) >= tolM12 * abs(Data(J)-Data(I))
                    % ... Schwingspiel gefunden
                    counter = counter + 1;
                    pz = pz + 1;
                    % ... Test ausgeschnittenes SSP
%                     [SUBDATA] = CutOutHyst(Data,I,J,K,NF);
%                     figure,grid on,plot(SUBDATA);
                    % ... Schwingspiele und Hysterese Schlie�index
                    AMP(1,pz) = counter;
                    AMP(2,pz) = K;
                    % ... Amplituden
                    AMP(3,pz) = 0.5* abs(Data(K)-Data(J));
                    AMP(4,pz) = 0.5* (Data(K)+Data(J));
                    AMP(5,pz) = 0.5* abs(Data(I)-Data(J));
                    AMP(6,pz) = 0.5* (Data(I)+Data(J));
                    IZ = IZ - 2;
                    % ... auschneiden Hysterese aus Zeitreihe
                    % durch manipulieren der Nachfolge und Vorg�nger Zeiger
                    % 1. Ende der Hysterese auf Dateiende 
%                     NF(VG(K)) = NF(1);
%                     VG(NF(1)) = VG(K);
%                     % 2. altes Datenende auf Anfang der Hysterese
%                     NF(1) = NF(I);
%                     VG(NF(I)) = 1;
                    % 3. Hystrese ausschneiden/�berspringen
                    NF(I) = K;
                    VG(K) = I;
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
                end % Ende Verzweigung �berpr�fung der Schwingweiten
            end % Ende Verzweigung UKP
            
        else
            % ... IZ <= IR (wird hier zsm behandelt)
            % ... Einlesen Wert aus Residuum
            J = RES(IZ);
            weiter = 0;
            % ... Pr�fe UKP
            if (Data(K)-Data(J))*Data(J) < 0
                % ... UKP
                % ... Pr�fe Memory 3
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












