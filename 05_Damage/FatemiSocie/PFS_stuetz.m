classdef PFS_stuetz < handle
% -------------------------------------------------------------------------
% Klassendefinition: Schädigungsparameter von Fatemi & Socie mit
% Wöhlerliniendefinition aus Stützstellen
% Zusammenfassung aller Definitionen/ Optionen & Funktionen für
% Schädigungsparameter PFS
%
%
%__________________________________________________________________________
% KONSTRUKTOR:
%
%    P = PFS(E,nu,kfs,Pfs_WS_stuetz,d1,d2,fak,varargin);
%      E,nu              - Elastizität
%      kfs               - Fatemi/Socie Normalspannungsempfindlichkeit 
%                    !!!!!!!!! geteilt durch Fließspannung
%     Pfs_WS_stuetz      - Stützstellen Pram-Wöhlerlinie bei definierter
%                          Lebensdauer für Werkstoff
%                          Pfs_stuetz bei N = 1000
%         d1,d2          - Neigungen PWL
%                          d1 für N < 1000
%                          d2 für N > 1000
%          fak           - Faktor Pfs_bauteil = 1/fak * Pram_WS_stuetz
%                          Zusammenfassung von Sicherheits & Stützfaktoren
%      varargin          - variabler Input
%                       'VarName1',var1,'VarName2',var2,...)
%
% mit varargin können alle anderen zusätzlichen public Eigenschaften
% gesetzt werden
%__________________________________________________________________________
%
% EIGENSCHAFTEN:
% Public:
%      cc          - counting channel, Zeiger auf die Zeile in Spannungen
%                    und dehnungen die mit HCM gezählt wird
%                            cc = 12 gamma_xz (Mode II an Punkt A)
%                    default cc = 10 gamma_xy (Mode III an Punkt A)
%  N_stuetz        - Stützstelle Lebensdauer
%                    default N_stuetz = 1000
%  ND_stuetz       - Stützstelle Lebensdauer Dauerfestigkeit
%                    default N_stuetz = 1e30
%  d1,d2           - Neigungen PWL
%                          d1 für P > P_stuetz
%    fak           - Faktor Pram_bauteil = 1/fak * Pram_WS_stuetz
%                          Zusammenfassung von Sicherheits & Stützfaktoren
% kfs              - Fatemi/Socie Normalspannungsempfindlichkeit 
%                    !!!!!!!!! geteilt durch Fließspannung
% sigF             - Fließspannung, mögliche Definition sigF = (Rm+Rp02)/2
%   ND             - Dauerfestigkeit Basis WL
%   E              - E Modul
%  nu              - Querdehnzahl
% miner            - option zum definieren der P-WL
%                    (0=elementar) (1=modifiziert) (sonst = orig)
%
% Protected:
%   PN             - Logaritmierte Werte der P-WL
%   PD             - Schädigungsparameter Dauerfestigkeit
% dauerfest        - Dummy Wert für dauerfestigkeit
% G                - Schubmodul (aus E & nu bestimmt)
% Name             - Name des Parameters (nur zum generieren von Output)
% Pfs_Bau_stuetz  - Stuetzstelle Bauteil
% Pfs_BauD_stuetz - Stuetzstelle Bauteil Dauerfest
% d3               - Steigung nach der Dauerfestigkeit
%__________________________________________________________________________
% FUNKTIONEN:
% Normal:
%   pfs             -> berechnet Schädigungsparameter
%   rainflow        -> rainflowzählung mit hcm
%   hcm             -> hcm zählung 
%   lebendauer      -> berechnet Lebensdauer
%   damage_akk      -> Schadensakkumulation
% 
% Static:
%__________________________________________________________________________
% EXTRERNE FUNKTIONEN:
%__________________________________________________________________________
% -------------------------------------------------------------------------
% AUTHOR: JAN KRAFT 
% STAND: Februar 2021
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% EIGENSCHAFTEN
% -------------------------------------------------------------------------
% PUBLIC:
   properties % Notwendige Größen für Kerbnäherung
      % Primäre Schädigungsvariable
      cc {mustBeNumeric} = 10;
      % statisch
      E {mustBeNumeric} = NaN;
      nu {mustBeNumeric} = NaN;
      sigF {mustBeNumeric} = NaN;
      % PWL
      Pfs_WS_stuetz {mustBeNumeric} = NaN;
      Pfs_WSD_stuetz {mustBeNumeric} = 1e-20; % default miner elementar
      d1 {mustBeNumeric} = NaN;
      d2 {mustBeNumeric} = NaN;
      fak   {mustBeNumeric} = 1;    
      N_stuetz {mustBeNumeric} = 1000;   
      ND_stuetz = NaN;
      % Minerregel
      miner {mustBeNumeric} = 0;
      % Fatemi Socie konstante
      kfs {mustBeNumeric} = NaN;
   end % Ende EIGENSCHAFTEN Public 
   
% PROTECTED
   properties (SetAccess = protected) % Größen werden automatisch gesetzt
      PN {mustBeNumeric} = NaN;
      PD {mustBeNumeric} = NaN; % Schädigungsparameter Dauerfestigkeit
      dauerfest = 1e20;         % definition für dauerfest
      G = NaN;                  % Schubmodul      
      Name = 'PFS_st';             % Name des Parameters
      Pfs_Bau_stuetz {mustBeNumeric} = NaN; % Stuetzstelle Bauteil
      Pfs_BauD_stuetz {mustBeNumeric} = NaN; % Stuetzstelle Bauteil
      d3 {mustBeNumeric} = NaN;              % Steigung nach Dauerfestigkeit
   end % Ende EIGENSCHAFTEN protected 
   
% -------------------------------------------------------------------------
% KONSTRUKTOR
% ------------------------------------------------------------------------- 
methods
    
    function obj = PFS_stuetz(E,nu,kfs,Pfs_WS_stuetz,d1,d2,fak,varargin)
        %------------------------------------------------------------------
        % Konstruktor der Klasse
        % INPUT:
        %      E,nu              - Elastizität
        %      sigF              - Fließspannung
        %     Pfs_WS_stuetz      - Stützstellen Pram-Wöhlerlinie bei definierter
        %                          Lebensdauer für Werkstoff
        %                          Pram_stuetz bei N = 1000
        %         d1,d2          - Neigungen PWL
        %                          d1 für N < 1000
        %                          d2 für N > 1000
        %          fak           - Faktor Pram_bauteil = fak * Pram_WS_stuetz
        %                          Zusammenfassung von Sicherheits & Stützfaktoren
        %           varargin     - variabler Input
        %                       'VarName1',var1,'VarName2',var2,...)
        % OUTPUT:
        %    obj - objekt der Klasse
        %------------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % Setze Felder
        % Elastisch
        obj.E = E;
        obj.nu = nu;
        obj.G = E/(2*(1+nu));
        obj.kfs = kfs;
        % DWL
        obj.Pfs_WS_stuetz = Pfs_WS_stuetz;
        obj.d1 = d1;
        obj.d2 = d2;
        obj.fak = fak;
     
        % -----------------------------------------------------------------
        % Variabler Input
        nvarargin = length(varargin);
        if nvarargin > 0
            for i = 2 : 2 : nvarargin
                varname = varargin{i-1};
                varvalue = varargin{i};
                if strcmp(varname,'miner')
                    obj.miner = varvalue;
                elseif strcmp(varname,'cc')
                    obj.cc = varvalue;
                elseif strcmp(varname,'kfs')
                    obj.kfs = varvalue;
                elseif strcmp(varname,'ND_stuetz')
                    obj.ND_stuetz = varvalue;
                elseif strcmp(varname,'N_stuetz')
                    obj.N_stuetz = varvalue;
                elseif strcmp(varname,'Pfs_WSD_stuetz')
                    obj.Pfs_WSD_stuetz = varvalue;
                else
                    msg = ['Die Eigenschaft ',varname,' wurde nicht erkannt'];
                    warning(msg);
                end
            end
        end
        
        
        % -----------------------------------------------------------------
        % bestimme P-WL
        % ... Stützstellen Bauteil
        obj.Pfs_Bau_stuetz = 1/fak * obj.Pfs_WS_stuetz;
        obj.Pfs_BauD_stuetz = 1/fak * obj.Pfs_WSD_stuetz;
        % ... Steigung nach der Dauerfestigkeit
        if obj.miner == 0 % elementar
            obj.d3 = d2;
            obj.PD = 1e-20;
        elseif obj.miner == 1 % modifiziert
            obj.d3 = d2/(2-d2);
            obj.PD = 1e-20;
        else % original
            obj.d3 = -1e-20;
            obj.PD = obj.Pfs_BauD_stuetz;
        end
        % ... Anpassen Dauerfestigkeit
        obj.ND_stuetz = obj.N_stuetz * (obj.Pfs_BauD_stuetz/obj.Pfs_Bau_stuetz)^(1/obj.d2);
        
       end % Ende Konstrutor
       
    end % Ende KONSTRUKTOREN Public   
   
% -------------------------------------------------------------------------
% METHODEN
% ------------------------------------------------------------------------- 
    methods
        % ... schädigungsparameter
        function P = pfs(obj,DATA,I,J,K)
            % Funktion zum berechnen des Schädigungsparameters von Fatemi/Socie
            %
            % INPUT:
            %  DATA      -> Spannungen und Dehnungen und Druchläufe
            %   I        -> Index Start Hysterese
            %   J        -> Index UKP Hysterese
            %   K        -> Index Ende Hysterese
            %
            % OUTPUT:
            %  P         -> Schädigungsparameter des gezählten Ssp. nach Fatemi/Socie
            %
            %__________________________________________________________________________
            
            % -------------------------------------------------------------------------
            % Indizes der Hysterese festlegen (unterscheide Auf- und Absteigende Äste)
            % und stehender/hängender Hysteresen
            % if DATA(cc,I) < DATA(cc,J) % Start der Hysterese kleine als UKP
            %     % ... stehende Hysterese
            %     iu_auf = I;
            %     io_auf = J;
            %     io_ab = J;
            %     iu_ab = K;
            % else
            %     % ... hängende Hysterese
            %     io_ab = I;
            %     iu_ab = J;
            %     iu_auf = J;
            %     io_auf = K;
            % end
            
            % -------------------------------------------------------------------------
            % Amplitude der Primären Schädiungsgröße (im normalfall gleitungen)
            % dgam = DATA(cc,io_auf) - DATA(cc,iu_auf);     % am aufsteigenden Ast
            % dgam = DATA(cc,io_ab) - DATA(cc,iu_ab);       % am absteigenden Ast
            gam_a = 0.5 * abs( DATA(obj.cc,I) - DATA(obj.cc,J) );        % willkürlich definiert
              
            % -------------------------------------------------------------------------
            % Maximale und minimale Normalspannung
            snmax = max(DATA(1,I:K));
            % imax = imax + I - 1;
            
            % -------------------------------------------------------------------------
            % Schädigungsparameter 
            P = gam_a * ( 1 + obj.kfs * snmax );
            
        end % Ende Berechnung Schädigungsparameter
        
        % ... rainflowzählung
        function P = rainflow(obj,sigepsfile,ndata,phi,psi)
            %==============================================================
            % rainflow zählung und Schädigungsbewertung
            %
            %==============================================================
            % Eingabeparameter:
            %
            % sigepsfile        - Verweis auf Datei mit diskreten
            %                     Lastpunkten für (elast-plast) Spannungen  
            %                     und Dehungen im Kerbkoordinantensystem
            %                     Data = [DLZ,S11,S22,S12,E11,E22,E33,G12]  
            %                     Data = [S,E,DLZ];
            %                     (double array e R^(13 x numdata))
            %                     DLZ -> Durchlaufzähler
            % ndata             - Anzahl Zeitpunkte in Lastfolge
            % phi,psi           - Winkel kritische Ebene (in Rad)
            %
            %==============================================================
            % Rückgabewerte:
            %
            % P                 - Schädigungsparameter der gezählten 
            %                     Schwingspiele
            %                     1. Zeile Schwingspielzähler
            %                     2. Zeile Schädigungsparameter
            %                     3. Zeile Durchlaufzähler
            %
            %==============================================================
            
            %==============================================================
            % Initialisierung
            % ... Ändere Zählvariable für 45°
            if psi == 45*pi/180
                obj.cc = 12;
            else
                obj.cc = 10;
            end
            % ... ein paar Definitionen
            nBuf = 10000;                                                  % Buffergröße
            nBufMax = 100e6;                                                % Buffergröße die nicht überschritten werden soll
            nSp = 13;                                                      % Anzahl der Spalten (1DLZ 3Spannungen 4 Dehnungen)
            datacounter = 0;                                               % Zählt wieviele Daten gelesen wurden
            % ... Öfnne Lastdatei & Lese Header
            fidLoad = fopen(sigepsfile,'r');                               % FileID
            % ... Vorbereiten Rainflow
            nRes = min(ndata,5e5);                                         % Beschränke Größe des Residuums
            nDam = min(ceil(ndata/2),5e5);                                 % Beschränke Größe der Outputvariable
            IR = 1;                                                        % Zeiger auf den letzten nicht schließfähigen Wert
            IZ = 1;                                                        % Zeiger auf den letzten noch schließfähigen Wert
            counter = 0;                                                   % Schleifenzähler
            RES = zeros(1,nRes); RES(1) = 1;                               % Residuum            
            Buf = zeros(nSp,nBuf);                                         % Speicher
            nextDat = 1;                                                   % Zeiger auf nächste freie Stelle im Buffer
            empty = nBuf;                                                  % Anzahl freie Stellen im Buffer
            Nf = 2:nBuf+1; Nf(nBuf) = 1;                                   % Zeiger auf Nachfolger im Buffer
            Vg = 0:nBuf-1; Vg(1) = nBuf;                                   % Zeiger auf Vorgänder im Buffer
            % ... Init Zeiger auf freie stelle in Schädigung
            P = zeros(3,nDam);                                             % Speicher für Schädigungsparameter
            p = 0;                                                         % zeiger auf nächste frei Spalte in P
            % ... Rotationsmatricen von Spannungen und Dehnungen in
            % kritische Ebene zu Drehen
            PS = transformstress(phi,psi);                                 % Dreht Spannungen
            PE = transformstrain(phi,psi);                                 % Dreht Dehnungen
            
            %==============================================================
            % Rainflow Zählung
            % Äußere Schleife über alle Werte in sigepsfile 
            status = 1;                                                    % Abbruchbedingung
            while status && empty > 1
                idx = NaN(1,nBuf);                                         % init Indices in Buffer in die beschrieben werden
                ndatared = 0;                                              % Anzahl (in diesem Durchlauf) gelesener Daten
                % 1. innere Schleife zum füllen des Buffers
                while status && Nf(nextDat) ~= RES(IR) 
                    % merke Index
                    idx(ndatared+1) = nextDat;
                    % Zeiger auf nächste freie Stelle Setzten
                    if datacounter < ndata
                        nextDat = Nf(nextDat);
                    else
                        status = 0;
                    end
                    % Inkrementiere Zähler
                    ndatared = ndatared + 1;
                    datacounter = datacounter + 1;                       
                end % Ende 1. innere Schleife
                idx(isnan(idx)) = [];                                      % Lösche dummy Werte
                % Lese Werte
                A = fread(fidLoad,[8,ndatared],'double');
                ncols = size(A,2);
                if isempty(A)
                    status = 0;
                else
                    % Drehe Spannungen
                    sig = PS * A(2:4,:);
                    % Drehe Dehnungen
                    eps = PE * A(5:8,:);
                    % Fülle Buffer                   
                    Buf(1:6,idx(1:ncols)) = sig;                           % Spannungen
                    Buf(7:12,idx(1:ncols)) = eps;                          % Dehnungen
                    Buf(13,idx(1:ncols)) = A(1,:);                         % Durchlaufzähler
                end
                % Buffer ist voll jetzt Zählen
                % 2. innere Schleife über die neuen Werte
                K = Nf(RES(IZ));                                           % Zeiger K auf Startindex setzten
                while K ~= nextDat

                    % Rainflow
                    [IZ,IR,RES,counter,Nf,Vg,P,p] = obj.hcm(...
                        K,Buf,Nf,Vg,IZ,IR,RES,counter,nextDat,P,p);
                    % Zeiger auf nächsten Wert setzten
                    K = Nf(K);
                end % Ende 2. innere Schleife
                
                % Freigewordene Stellen (geschlossene Hyst ermitteln)
                empty = 1;
                i = nextDat;
                while Nf(i) ~= RES(IR)
                    empty = empty + 1;
                    i = Nf(i);
                end
%                 fprintf('empty:%8i nBuf: %8i Nf(nextDat):%8i RES(IR):%8i IR:%i\n',empty,nBuf,Nf(nextDat),RES(IR),IR);
                % Falls der Buffer voll ist aber Daten nicht zu ende sind
                % überlauf anzeigen
                if empty <= 1 && 2*nBuf < nBufMax
%                     msg = 'Bufferüberlauf in Rainflow PFS:Buffergröße verdoppeln';
%                     warning(msg)
                    % verdopple Buffergröße
                    Buf = [Buf,zeros(nSp,nBuf)];                           % Neuer Buffer
                    Vg = [Vg,nBuf:2*nBuf-1];                               % Neuer Vorgänger
                    Vg(nBuf + 1) = nextDat;
                    Vg(Nf(nextDat)) = 2*nBuf; 
                    Nf = [Nf,nBuf+2:2*nBuf+1];                             % Neuer Nachfolger 
                    Nf(2*nBuf) = Nf(nextDat); 
                    Nf(nextDat) = nBuf + 1;         
                    empty = empty + nBuf;
                    nBuf = 2*nBuf;
                elseif empty <= 1
                    msg = 'Bufferüberlauf Programm Rainflow PFS\n';
                    error(msg);                          
                end
            end % Ende Schleife über alle Werte
            
            P = P(:,P(1,:) ~= 0);  % Aussortieren leere P Werte
            
            %==============================================================
            % Schließe Lastdatei
            fclose(fidLoad);
            
        end % Ende Rainflowzählung
        
        % ... hcm Algorithmus
        function [IZ,IR,RES,counter,Nf,Vg,P,p] = hcm(obj,...
                  K,Data,Nf,Vg,IZ,IR,RES,counter,nextDat,P,p)
            % -------------------------------------------------------------
            % Hilfsfunktion zum zyklenzählen
            % HCM Hysteresis Counting Method
            %
            % Der vorliegende Code wurde auf Grundlage von
            % RAINFLOW-HCM Ein Hyseresisschleifen-Zählalgorithmus auf
            % Werkstoffmechanischer Grundlage von U.H. Chlormann und 
            % T. Seeger aus dem Jahr 1985 implementiert
            % -------------------------------------------------------------
            % Toleranzen
            tolM12 = 0.99;                             % Toleranz für das Erkennen vom Memory 1 und 2
            tolM3 = 1.01;                              % Toleranz für das Erkennen vom Memory 3
            
            % Abbruchbedingung
            weiter = 1;
            
            % =============================================================
            % Rainflow
            while weiter
                % 2
                if IZ > IR % Vergleich der Zeiger
                    
                    % ... letzte Werte aus Residuum lesen
                    I = RES(IZ-1);
                    J = RES(IZ);
                    
                    % ... Prüfe ob letzter Wert UKP ist
                    if (Data(obj.cc,K)-Data(obj.cc,J))*(Data(obj.cc,J)-Data(obj.cc,I)) >= 0%-1e-20
                        % ... kein UKP
                        IZ = IZ - 1;
                        weiter = 1;
                        %                 !! GOTO 2
                    else
                        % ... UKP
                        % ... Prüfe Schwingweite größer als die letzte
                        if abs(Data(obj.cc,K)-Data(obj.cc,J)) >= tolM12 * abs(Data(obj.cc,J)-Data(obj.cc,I))
                            % ... Schwingspiel gefunden
                            counter = counter + 1;
                            % ... Schneide Schwingspiel aus
                            [SUBDATA,In,Jn,Kn] = CutOutHyst2(Data,I,J,K,Nf);
                            % ... Schädigungsrechnung
                            Pfs = obj.pfs(SUBDATA,In,Jn,Kn);                          
                            % ... Speichern Schwingspiel- & Durchlaufzähler
                            % wenn größer als Dauerfestigkeit
                            if Pfs >= obj.PD
                                p = p + 1;
                                P(:,p) = [counter;Pfs;Data(13,K)];
                            end
                            % ... Dekrementieren Zeiger
                            IZ = IZ - 2;
                            % ... auschneiden Hysterese aus Zeitreihe
                            % Ende der Hysterese auf altes Datenende
                            Nf(Vg(K)) = Nf(nextDat);
                            Vg(Nf(nextDat)) = Vg(K);
                            % Altes Datenende auf Anfang der Hysterese
                            Nf(nextDat) = Nf(I);
                            Vg(Nf(I)) = nextDat;
                            % Hysterese ausschneiden/überspringen
                            Nf(I) = K;
                            Vg(K) = I;
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
                    if (Data(obj.cc,K)-Data(obj.cc,J))*Data(obj.cc,J) < 0
                        % ... UKP
                        % ... Prüfe Memory 3
                        if abs(Data(obj.cc,K)) > tolM3 * abs(Data(obj.cc,J))
                            % ... Memory 3
                            IR = IR + 1;
                        end
                        IZ = IZ + 1;
                    end
                end % Ende Verzweigung Zeigervergleich
            end % Ende while Schleife
            
            RES(IZ) = K;
        end % Ende hcm
        
        % ... Lebensdauer Rechnung 
        function [DL,SSP,PDam] = lebensdauer(obj,P,ndl)
            % Funktion rechnet Lebensdauern aus Schädigungsparametern
            %
            % INPUT:
            % P              - Schädigungsparameter
            %                  1.Zeile Schwingspiele prim. Schädigungsvariable
            %                  2.Zeile Schädigungsparameter
            %                  3.Zeile Durchläufe
            % ndl            - Durchläufe der Lastfolge, die Simuliert
            %                  wurden
            %
            % OUTPUT:
            % DL             - Durchläufe
            % SSP            - Schwingspiele
            % P              - Erweitert um 
            %                  4. Zeile aktueller Schädigungsbeitrag
            %                  5. Zeile akkumilierte Schädigung
            %______________________________________________________________
            
            %--------------------------------------------------------------
            % Schädigungsrechnunge
            nP = size(P,2);            % Anzahl Werte in P 
            PDam = zeros(2,nP);        % Speicher für Schädigung
            Dsum = 0;                  % Schadenssumme
            Dlast = 0;                 % Schädigung letzter Durchlauf
%             ndl = ceil(max(P(3,:)));   % maximale Anzahl an Durchläufen
            idam = 0;                  % Zeiger auf den Wert an dem Dsum = 1
            ilast = 1;                 % Zeiger auf ersten Wert des letzten Durchlaufs
            for i = 1: nP
                % ... Schädigung aus WL
                [Dakt,Dsum] = obj.damage_akk(P(2,i),Dsum);
                PDam(:,i) = [Dakt;Dsum];
                % ... Schädigung des letzten Durchlaufsspeichern
                if P(3,i) > ndl - 1
                    Dlast = Dlast + Dakt;
                else
                    ilast = i;
                end
                % ... kaputt
                if Dsum >= 1
                    idam = i;
                    break;
                end
            end
            
                        
            % Anzahl Zyklen letzter Durchlauf
            if nP > 0
                cyclast = P(1,nP) - P(1,ilast);  % Schwingspiele im letzten DL
                nssp = P(1,nP);                  % Anzahl an Swingspielen
            else
                DL = obj.dauerfest;
                SSP = obj.dauerfest;
                return;
            end
            
            
            %--------------------------------------------------------------
            % Schauen obs kaputt is
            if Dsum >= 1
                % ... kaputt
                DL = P(3,idam);
                SSP = P(1,idam);
            else
                % ... nicht kaputt
                if Dlast > 0
                    % ... 1 = D + x*Dlast, x = noch zu ertragende
                    %                          Durchläufe
                    x = (1-Dsum)/Dlast;
%                     DL = P(3,nP) + x;
                    DL = ndl + x;
                    SSP = nssp + int8(x*cyclast); % rundet auf
                else
                    % ... Dauerfest
                    DL = obj.dauerfest;
                    SSP = obj.dauerfest;
                end
            end

        end % Ende Lebensdauerrechnung
        
        % ... Schadensakkumulation
        function [Dakt,Dsum] = damage_akk(obj,P,Dsum)
            % Funktion berechnet Schädigung aus Schädigungsparameter und
            % logaritmierter P-WL
            %
            % INPUT:
            % P          - (double) Schädigungsparameter
            % Dsum       - Schadenssumme
            %
            % OUTPUT:
            % Dakt       - aktuelle Schädigung aus P
            % Dsum       - aktualisierte Schadenssumme
            %______________________________________________________________
            % Init aktueller Schadensbeitrag
            Dakt = 0;

            % Sicherstellen P>0
            if P > 0
                % P > Ps ?
                if P > obj.Pfs_Bau_stuetz % Steigung 1/d1
                    Dakt = 1/obj.N_stuetz * (obj.Pfs_Bau_stuetz/P)^(1/obj.d1);
                elseif  P > obj.Pfs_BauD_stuetz % Steigung 1/d2
                    Dakt = 1/obj.N_stuetz * (obj.Pfs_Bau_stuetz/P)^(1/obj.d2);
                else % Unterhalb der Dauerfestigkeit Steigung 1/d3
                    Dakt = 1/obj.ND_stuetz * (obj.Pfs_BauD_stuetz/P)^(1/obj.d3);
                    % Abfangen NaN Werte die durch unendliche Lebensdauer
                    % enstehen kann
                    if isnan(Dakt)
                        Dakt = 0;
                    end
                end
            end

            % Akkumulation
            Dsum = Dsum + Dakt;
        end % Ende Schadensakkumulation


        % ... plote P-WL
        function plotPWL(obj,varargin)
            % Plote Wöhlerlinie
            % varargin = ax
            if nargin == 1
                figure, grid on, hold on
                ax = gca;
            else
                ax = varargin{1};
            end
            plot([1, obj.N_stuetz, obj.ND_stuetz, 1e20], ...
            [(1/obj.N_stuetz)^obj.d1*obj.Pfs_Bau_stuetz, obj.Pfs_Bau_stuetz, obj.Pfs_BauD_stuetz,(1e20/obj.ND_stuetz)^obj.d3*obj.Pfs_BauD_stuetz])
            set(ax,'XScale','log','YScale','log');
            set(ax,'XLim',[1 1e7])
            xlabel('N'), ylabel('P')
            grid on
        end  % Ende plot PWL


    end % Ende Methoden

% -------------------------------------------------------------------------
% STATISCHE METHODEN
% ------------------------------------------------------------------------- 
    methods (Static)
    end % Ende statische Methoden

end % Ende Klassendefinition PRAM
