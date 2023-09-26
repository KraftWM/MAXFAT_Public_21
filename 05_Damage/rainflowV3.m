function P = rainflowV3(Data,DMGs,psi)
% HCM Hysteresis Counting Method um Amplituden in einem Signal zu finden
% 
% Im Residuum werden keine Werte sondern nur Spaltenindices
%       gespeichert 
%
% Gezählt werden die drei Schnittdehnungen
%  exx -> Mode I     (cc = 7 )
%  gxz -> Mode II    (cc = 12)
%  gxy -> Mode III   (cc = 10)
%
% !!!
% Es werden nur die Umkehrpunkte gezählt (wegen Geschwindigkeit). Dadurch
% können die tatsächlichen Schließzeitpunkte im Vorfeld nicht bestimmt
% werden. Es kann so einer geringfügigen Veränderung der Lebensdauer im
% Kurzrissmodell kommen, weil Schwingspiele evt. doch früher schließen. Der
% Einfluss ist allerdings sehr gering bis nicht existent.
% !!!
%
% Der vorliegende Code wurde nach
% RAINFLOW-HCM Ein Hyseresisschleifen-Zählalgorithmus auf
% Werkstoffmechanischer Grundlage von U.H. Chlormann und T. Seeger aus dem
% Jahr 1985 implementiert
%
%==========================================================================
% INPUT:
%
%
% Aus Kerbsimulation
% Data              - Diskrete Lastpunkt für Spannungen/Dehungen/Lasten im
%                     Koordinatensystem der kritischen Ebene
%                     Data(1:6,:)  -> SXX,SYY,SZZ,SXY,SYZ,SXZ
%                     Data(7:12,:) -> EXX,EYY,EZZ,GXY,GYZ,GXZ
%                     Data(13,:)   -> DLZ
%
% DMGs              - cell array mit Objekten einer 
%                     Schädigungsparameterklasse, aktuell kompatibel mit:
%                     PRAM,PRAM_stuetz
%                     PFS,PFS_stuetz
%                     PZ,PZ2,PZ2_stuetz
%
%==========================================================================

%
% OUTPUT:
% P                 - cell array Schädigungsparametern, P{i} entspricht 
%                     dem Input in die Funktion DGMs{i}.lebensdauer
%==========================================================================
% Erstellt von Jan Kraft
% Version 2.0 September 2021
%==========================================================================


%% Verwalte Schädigungsmodelle
ndata = size(Data,2);
nDMGs = size(DMGs,2);
p = zeros(1,nDMGs);
P = cell(1,nDMGs);
% Verwalte DMGs (vllt muss nicht alles gezählt werden
mode1 = zeros(1,nDMGs);    % Mode 1 (exx) nicht zählen 
mode2 = zeros(1,nDMGs);    % Mode 2 (gxz) nicht zählen 
mode3 = zeros(1,nDMGs);    % Mode 3 (gxy) nicht zählen 
for i = 1:length(DMGs)
    if any(strcmp(DMGs{i}.Name,{'PZ','PZ_st','PZ_neu'}))
        mode1(i) = 1;
        mode2(i) = 1;
        mode3(i) = 1;
        P{i} = zeros(4, min(ndata,5e5));
    elseif any(strcmp(DMGs{i}.Name,{'PFS_mod'}))
        % Für Fatemi/Socie in Ebenen psi~=0 wird nicht mit Rainflow sondern
        % mit Modified Wang Brown gezählt
        % Speicher für P{i} wird tortzdem initialisiert
        % -> kein Mode wird Eingeschaltet
        if abs(psi) < 1e-10
            if DMGs{i}.cc == 12
                mode2(i) = 1;
            elseif DMGs{i}.cc == 10
                mode3(i) = 1;
            end
        end
        P{i} = zeros(4, min(ndata,5e5));
    elseif any(strcmp(DMGs{i}.Name,{'PFS_st','PFS','PFS_var'}))
        % Umschalten Fatemi/Socie in psi = 45° auf Mode II
        if psi == 45*pi/180
           DMGs{i}.cc = 12;
        else
            DMGs{i}.cc = 10;
        end
        if DMGs{i}.cc == 12
            mode2(i) = 1;
        elseif DMGs{i}.cc == 10
            mode3(i) = 1;
        end
        P{i} = zeros(3, min(ndata,5e5));
    elseif any(strcmp(DMGs{i}.Name,{'PRAM','PRAM_st'}))
        if DMGs{i}.cc == 7
            mode1(i) = 1;
        end
        P{i} = zeros(3, min(ndata,5e5));
    end
end

%% Ermittle & Sortiere Umkehrpunkte
UKP = [];                                  % Init Speicher für Umkehrpunkte
nMode1 = 0;
nMode2 = 0;
nMode3 = 0;
% ... Mode 1
if any(mode1)
    % ... filter UKP
    [~,UKP1] = filterUKP(Data(7,:));
    nMode1 = size(UKP1,1);
    UKP = [UKP1,repelem(1,nMode1,1)];
end
% ... Mode 2
if any(mode2)
    % ... filter UKP
    [~,UKP2] = filterUKP(Data(12,:));
    nMode2 = size(UKP2,1);
    UKP = [UKP;[UKP2,repelem(2,nMode2,1)]];
end  
% ... Mode 3
if any(mode3)
    % ... filter UKP
    [~,UKP3] = filterUKP(Data(10,:));
    nMode3 = size(UKP3,1);
    UKP = [UKP;[UKP3,repelem(3,nMode3,1)]];
end
% ... sortieren
UKP = sortrows(UKP);
nUKP = size(UKP,1);


%% Initialisiere Rainflow Variablen
IR1 = 1;                                                                   % Zeiger auf den letzten nicht schließfähigen Wert (Mode I)
IR2 = 1;                                                                   % (Mode II)
IR3 = 1;                                                                   % (Mode III)
IZ1 = 1;                                                                   % Zeiger auf den letzten noch schließfähigen Wert (Mode I)
IZ2 = 1;                                                                   % (Mode II)
IZ3 = 1;                                                                   % (Mode III)
counter1 = 0;                                                              % Schleifenzähler (Mode I)
counter2 = 0;                                                              % (Mode II)
counter3 = 0;                                                              % (Mode III)
RES1 = zeros(1,nMode1); RES1(1) = 1;                                       % Residuum (Mode I)
RES2 = zeros(1,nMode2); RES2(1) = 1;                                       % (Mode II)
RES3 = zeros(1,nMode3); RES3(1) = 1;                                       % (Mode III)
keep1 = true(1,ndata);                                                     % Bool zum auschneiden von Hysteresen (Mode I)
keep2 = true(1,ndata);                                                     % (Mode II)
keep3 = true(1,ndata);                                                     % (Mode III)

%% Rainflow
for k = 1 : nUKP
    
    % ... Aktueller Mode
    aktmode = UKP(k,2);
    
    % ... Mode 1
    if aktmode == 1 && any(mode1)
      [IZ1,IR1,RES1,counter1,keep1,P,p] = hcm(...
            UKP(k,1),Data,IZ1,IR1,RES1,counter1,7,1,keep1,P,p,DMGs);  

    
    % ... Mode 2
    elseif aktmode == 2 && any(mode2)
      [IZ2,IR2,RES2,counter2,keep2,P,p] = hcm(...
            UKP(k,1),Data,IZ2,IR2,RES2,counter2,12,2,keep2,P,p,DMGs);  
    
    
    % ... Mode 3
    elseif aktmode == 3 && any(mode3)
      [IZ3,IR3,RES3,counter3,keep3,P,p] = hcm(...
            UKP(k,1),Data,IZ3,IR3,RES3,counter3,10,3,keep3,P,p,DMGs);  
    end
    
    
end % End Schleife über alle Umkehrpunkte

%% Aufräumen
% clear Memory für Dehnungsvorgeschichte im Kurzrissmodell & aussortieren
% der leeren P-Werte
for i = 1:length(DMGs)
    if any(strcmp(DMGs{i}.Name,{'PZ','PZ_st','PZ_neu'}))
        DMGs{i}.exmax = 0;
        DMGs{i}.exmin = 0;
        DMGs{i}.exop_alt = 0;
        DMGs{i}.exop_ein = 0;
        P{i} = sortrows(P{i}(:,1:p(i))',4)';
    else
        P{i} = P{i}(:,1:p(i));
    end
end

end % Ende Hauptfunktion

% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
%                                                                         %
%                   Hilfsfunktion HCM für einzelne Werte                  %
%                                                                         %
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %

function [IZ,IR,RES,counter,keep,P,p] = hcm(...
    K,Data,IZ,IR,RES,counter,cc,mode,keep,P,p,DMGs)
% -----------------------------------------------------------------
% Hilfsfunktion zum zyklenzählen
% HCM Hysteresis Counting Method
%
% Der vorliegende Code wurde auf Grundlage von
% RAINFLOW-HCM Ein Hyseresisschleifen-Zählalgorithmus auf
% Werkstoffmechanischer Grundlage von U.H. Chlormann und
% T. Seeger aus dem Jahr 1985 implementiert
% -----------------------------------------------------------------

% Toleranzen
tolM12 = 0.99;                             % Toleranz für das Erkennen vom Memory 1 und 2
tolM3 = 1.01;                              % Toleranz für das Erkennen vom Memory 3

% Abbruchbedingung
weiter = 1;

% ========================================================================
% Rainflow
while weiter
    % 2
    if IZ > IR % Vergleich der Zeiger
        
        % ... letzte Werte aus Residuum lesen
        I = RES(IZ-1);
        J = RES(IZ);
        
        % ... Prüfe ob letzter Wert UKP ist
        if (Data(cc,K)-Data(cc,J))*(Data(cc,J)-Data(cc,I)) >= 0%-1e-20
            % ... kein UKP
            IZ = IZ - 1;
            weiter = 1;
            %                 !! GOTO 2
        else
            % ... UKP
            % ... Prüfe Schwingweite größer als die letzte
            if abs(Data(cc,K)-Data(cc,J)) >= tolM12 * abs(Data(cc,J)-Data(cc,I))
                % ... Schwingspiel gefunden
                counter = counter + 1;
                % Finde echtes K
                DJI = tolM12 * abs(Data(cc,J)-Data(cc,I));
                Khat = K - 1;
                while Khat > J && abs(Data(cc,Khat)-Data(cc,J)) >= DJI
                    Khat = Khat - 1;
                end
                Khat = Khat + 1;
                SubData = Data(:,I:Khat);
                SubData = SubData(:,keep(I:Khat));
                In = 1;
                Kn = size(SubData,2);
                Jn = length(find(keep(I:J)));
                keep(I+1:Khat-1) = false;
              
                % ---------------------------------------------------------
                % Unterscheide einzelnen Schädigungsparameter   
                for i = 1:length(DMGs)
                    % ... Schädigungsrechnung
                    switch DMGs{i}.Name                        
                        % Kurzrissmodell
                        case {'PZ','PZ_st'}                           
                            % ... Unterscheide Moden
                            if mode == 1
                                Pz = DMGs{i}.kurzriss_mode1(SubData,In,Jn,Kn);             
                            elseif mode == 2 
                                Pz = DMGs{i}.kurzriss_mode23( mode,SubData,In,Jn,Kn);                    
                            elseif mode == 3                               
                                Pz = DMGs{i}.kurzriss_mode23( mode,SubData,In,Jn,Kn);                               
                            end
                            % ... merken durchlauf am letzten Punkt der Hyst
                            if Pz >=  DMGs{i}.PZD0 * 0.1
                                p(i) = p(i) + 1;
                                P{i}(:,p(i)) = [mode;counter;Pz;Data(13,Khat)];
                            end
                            % ... Rissöffnungsdehnung abklingen lassen !!!
                            DMGs{i}.exop_alt =  DMGs{i}.exop_ein + (DMGs{i}.exop_alt - DMGs{i}.exop_ein)...
                                * exp(-DMGs{i}.FA/DMGs{i}.Q*Pz^DMGs{i}.mJ);
                        % Kurzrissmodell neu
                        case {'PZ_neu'}                           
                            % ... Unterscheide Moden
                            if mode == 1
                                Pz = DMGs{i}.kurzriss_mode1(SubData,In,Jn,Kn);             
                            elseif mode == 2 
                                Pz = DMGs{i}.kurzriss_mode23( mode,SubData,In,Jn,Kn);                    
                            elseif mode == 3                               
                                Pz = DMGs{i}.kurzriss_mode23( mode,SubData,In,Jn,Kn);                               
                            end
                            % ... merken durchlauf am letzten Punkt der Hyst
                            if Pz >=  DMGs{i}.PZD0 * 0.1
                                p(i) = p(i) + 1;
                                P{i}(:,p(i)) = [mode;counter;Pz;Data(13,Khat)];
                            end
                            % ... Rissöffnungsdehnung abklingen lassen !!!
                            DMGs{i}.exop_alt =  DMGs{i}.exop_ein + (DMGs{i}.exop_alt - DMGs{i}.exop_ein)...
                                * exp(-DMGs{i}.FA/DMGs{i}.Q(mode)*Pz^DMGs{i}.mJ(mode));
%                             fprintf('%.5f\n',DMGs{i}.exop_alt)
                        % Fatemi/Socie
                        case {'PFS', 'PFS_var','PFS_st'}
                            % ... Nur bewerten wenn Mode richtig is
                            if cc == DMGs{i}.cc
                                % ... berechne Schädigung
                                Pfs = DMGs{i}.pfs(SubData,In,Jn,Kn);
                                % ... Merke Zähler, Pfs und Schließzeitpunkt
                                if Pfs >= DMGs{i}.PD
                                    p(i) = p(i) + 1;
                                    P{i}(:,p(i)) = [counter;Pfs;Data(13,Khat)];
                                end
                            end
                        % Fatemi/Socie mit MWB Methode
                        case {'PFS_mod'}
                            % ... Nur bewerten wenn Mode richtig is
                            if cc == DMGs{i}.cc
                                % ... berechne Schädigung
                                Pfs = DMGs{i}.pfs(SubData,In,Jn,Kn);
                                % ... Merke Zähler, Pfs und Schließzeitpunkt
                                if Pfs >= DMGs{i}.PD
                                    p(i) = p(i) + 1;
                                    P{i}(:,p(i)) = [counter;Pfs;Data(13,Khat);1];
                                end
                            end                                
                        % PRAM    
                        case {'PRAM','PRAM_st'}
                            % ... Nur bewerten wenn Mode richtig is
                            if cc == DMGs{i}.cc
                                % ... berechne Schädigung
                                Pram = DMGs{i}.pram(SubData,In,Jn,Kn);
                                % ... Merke Zähler, Pfs und Schließzeitpunkt
                                if Pram >= DMGs{i}.PD
                                    p(i) = p(i) + 1;
                                    P{i}(:,p(i)) = [counter;Pram;Data(13,Khat)];
                                end
                            end
                            
                        otherwise
                            msg = ['Fehler in Rainflow: ',...
                                    'angegebenen Schädigungsparameter ',...
                                    'nicht erkannt'];
                           error(msg)
                    end % Ende Fallunterscheidung Schädigungsparameter
                    
                end % Ende Schleife über Schädigungsparameter
                % ---------------------------------------------------------
                
                % ... Dekrementieren Zeiger
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
            end % Ende Verzweigung überprüfung der Schwingweiten
        end % Ende Verzweigung UKP
        
    else
        % ... IZ <= IR (wird hier zsm behandelt)
        % ... Einlesen Wert aus Residuum
        J = RES(IZ);
        weiter = 0;
        % ... Prüfe UKP
        if (Data(cc,K)-Data(cc,J))*Data(cc,J) < 0
            % ... UKP
            % ... Prüfe Memory 3
            if abs(Data(cc,K)) > tolM3 * abs(Data(cc,J))
                % ... Memory 3
                IR = IR + 1;
            end
            IZ = IZ + 1;
        end
    end % Ende Verzweigung Zeigervergleich
end % Ende while Schleife

RES(IZ) = K;
end % Ende hcm






