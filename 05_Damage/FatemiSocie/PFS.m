classdef PFS < handle
% -------------------------------------------------------------------------
% Klassendefinition: Schädigungsparameter von Fatemi & Socie
% Zusammenfassung aller Definitionen/ Optionen & Funktionen für
% Schädigungsparameter PFS
%
%
%__________________________________________________________________________
% KONSTRUKTOR:
%
%    P = PFS(E,nu,sigF,sf,ef,b,c,ND,varargin);
%      E,nu              - Elastizität
%      sigF              - Fließspannung
%    sf,ef,b,c,ND        - Parameter Dehnungswöhlerlinie
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
%                    und Dehnungen die mit HCM gezählt wird
%                            cc = 12 gamma_xz (Mode II an Punkt A)
%                    default cc = 10 gamma_xy (Mode III an Punkt A)
% pwlopt           - option zum bestimmen der P-Wöhlerlinie
%                    'dwl' -> aus Dehnungswöhlerlinie (default)
%                    'gwl' -> aus Gleitungswöhlerlinie
% kfsopt           - option zum behandeln des k parameters
%                    'konst' -> konstantes k (default)
%                    'var'   -> k als Funktion der Lebensdauer k=k(N)
% sf,ef,b,c        - Parameter Dehnungswöhlerlinie
% tf,gf,b0,c0      - Parameter Gleitungswöhlerlinie
% fwt              - Schubwechselfestigkeit zu Wechselfestigkeit
%                    default fwt = 1/sqrt(3) -> Mises
% kfs              - Fatemi/Socie Normalspannungsempfindlichkeit 
%                    !!!!!!!!! geteilt durch Fließspannung
% sigF             - Fließspannung, mögliche Definition sigF = (Rm+Rp02)/2
%   ND             - Dauerfestigkeit Basis WL
%   E              - E Modul
%  nu              - Querdehnzahl
% nst              - Stützziffer
% miner            - option zum definieren der P-WL
%                    (0=elementar) (1=modifiziert) (sonst = orig)
%
% Protected:
%   PN             - Logaritmierte Werte der P-WL
%   PD             - Schädigungsparameter Dauerfestigkeit
%   q              - Interpolationsparameter zwischen GEH und NH
% dauerfest        - Dummy Wert für dauerfestigkeit
% G                - Schubmodul (aus E & nu bestimmt)
% toleranz         - tolearanz für newtoniteration
% N0               - Startwert Newtoniteration
% Name             - Name des Parameters (nur zum generieren von Output)
%__________________________________________________________________________
% FUNKTIONEN:
% Normal:
%   pfs             -> berechnet Schädigungsparameter
%   rainflow        -> Rainflowzählung mit hcm
%   hcm             -> hcm Zählung 
%   lebendauer      -> berechnet Lebensdauer
%   damage_akk      -> Schadensakkumulation
% 
% Static:
%   dwl2pwl         -> berechne PFS-Wöhlerlinie aus Dehnungswöhlerlinie
%   gwl2pwl         -> berechne PFS-Wöhlerlinie aus Gleitungswöhlerline
% bestimmekfs       -> Konstanten k Parameter bestimmen 
% zielfun           -> Zielfunktion für Newtonverfahren für rechnung mit
%                      variablem k
% Diffzielfun       -> analytische Ableitung Zielfunktion
% skaliereGEHundNH  -> Interpolation zwischen GEH und NH 
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
      % optionen 
      pwlopt = 'dwl';
      kfsopt = 'konst';
      % Primäre Schädigungsvariable
      cc {mustBeNumeric} = 10;
      % statisch
      E {mustBeNumeric} = NaN;
      nu {mustBeNumeric} = NaN;
      sigF {mustBeNumeric} = NaN;
      % DWL
      sf {mustBeNumeric} = NaN;
      ef {mustBeNumeric} = NaN;
      b {mustBeNumeric} = NaN;
      c {mustBeNumeric} = NaN;
      ND {mustBeNumeric} = NaN;
      % GWL
      tf {mustBeNumeric} = NaN;
      gf {mustBeNumeric} = NaN;
      bg {mustBeNumeric} = NaN;
      cg {mustBeNumeric} = NaN;
      % Schubwechselfestigkeit zu Wechselfestigkeit
      fwt {mustBeNumeric} = 1/sqrt(3);
      % Minerregel
      miner {mustBeNumeric} = 0;
      % Stützziffer
      nst {mustBeNumeric} = 1;
      % Fatemi Socie konstante
      kfs {mustBeNumeric} = NaN;
   end % Ende EIGENSCHAFTEN Public 
   
% PROTECTED
   properties (SetAccess = protected) % Größen werden automatisch gesetzt
      PN {mustBeNumeric} = NaN;
      PD {mustBeNumeric} = NaN; % Schädigungsparameter Dauerfestigkeit
      q {mustBeNumeric} = NaN;  % Interpolationsparameter
      dauerfest = 1e20;         % definition für dauerfest
      G = NaN;                  % Schubmodul
      toleranz = 1e-12;          % Toleranz
      dN = 1;                   % Schrittweite zentrale Differenzen
      maxiter = 100;            % Maximale iterationen Newtonverfahren
      Name = 'PFS';             % Name des Parameters
   end % Ende EIGENSCHAFTEN protected 
   
% -------------------------------------------------------------------------
% KONSTRUKTOR
% ------------------------------------------------------------------------- 
methods
    
    function obj = PFS(E,nu,sigF,sf,ef,b,c,ND,varargin)
        %------------------------------------------------------------------
        % Konstruktor der Klasse
        % INPUT:
        %      E,nu,sigF         - Elastizität
        %    sf,ef,b,c,ND        - Parameter Dehnungswöhlerlinie
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
        obj.sigF = sigF;
        % DWL
        obj.sf = sf;
        obj.ef = ef;
        obj.b = b;
        obj.c = c;
        % Dauerfestigkeit
        obj.ND = ND;
     
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
                elseif strcmp(varname,'pwlopt')
                    obj.pwlopt = varvalue;
                elseif strcmp(varname,'kfsopt')
                    obj.kfsopt = varvalue;
                elseif strcmp(varname,'kfs')
                    obj.kfs = varvalue;
                elseif strcmp(varname,'nst')
                    obj.nst = varvalue;
                elseif strcmp(varname,'tf')
                    obj.tf = varvalue;
                elseif strcmp(varname,'gf')
                    obj.gf = varvalue;
                elseif strcmp(varname,'bg')
                    obj.bg = varvalue;
                elseif strcmp(varname,'cg')
                    obj.cg = varvalue;
                elseif strcmp(varname,'fwt')
                    obj.fwt = varvalue;
                else
                    msg = ['Die Eigenschaft ',varname,' wurde nicht erkannt'];
                    warning(msg);
                end
            end
        end
        
        % -----------------------------------------------------------------
        % Grenzen fwt
        if obj.fwt < 1/sqrt(3)
            obj.fwt = 1/sqrt(3);
        elseif obj.fwt > 1
            obj.fwt = 1;
        end
        
        % -----------------------------------------------------------------
        % Setze Berechnungsoptionen 
        % (abfangen nichtausführbarer Pfade, kfsopt = 'var' und pwlopt =
        % 'gwl' nur wenn Gleitungswöhlerlinie bekannt ist
        if any(isnan([obj.tf,obj.gf,obj.bg,obj.cg])) 
            if strcmp(obj.kfsopt,'var')
                msg = ['variables k nur bei bekannter GleitungsWL. ',...
                       'kfsopt = konst gesetzt'];
                warning(msg);
                obj.kfsopt = 'konst';
            end
            if strcmp(obj.pwlopt,'gwl')
                msg = ['Schädigung aus GleitungsWL nur bei bekannter ',...
                       'GleitungsWL. pwlopt = dwl gesetzt'];
                warning(msg)
                obj.pwlopt = 'dwl';
            end
        end
        
        % -----------------------------------------------------------------
        % Bestimme kfs falls nicht selbst gesetzt 
        if isnan(obj.kfs)
            % Mit abgeschätzter GWL
            if any(isnan([obj.tf,obj.gf,obj.bg,obj.cg])) 
                % Interpolationsparameter
                obj.q = obj.skaliereGEHundNH(obj.fwt);
                % Abschätzen GWL aus Normaldehnungshypothese
                [tf,gf,bg,cg] = GleitungsWL(obj.sf,obj.ef,obj.b,obj.c,3,obj.nu);
                [~,kfsNH] = obj.bestimmekfs(obj.sf,obj.ef,obj.b,obj.c,...
                                      tf,gf,bg,cg,...
                                      obj.sigF,obj.nu,obj.E);
                % Abschätzen GWL aus Gestaltänderungsenergiehypothese
                [tf,gf,bg,cg] = GleitungsWL(obj.sf,obj.ef,obj.b,obj.c,1);
                [~,kfsGEH] = obj.bestimmekfs(obj.sf,obj.ef,obj.b,obj.c,...
                                      tf,gf,bg,cg,...
                                      obj.sigF,obj.nu,obj.E);
                % Interpolation
                kfs = obj.q * kfsNH + (1-obj.q) * kfsGEH;
                obj.kfs = kfs/obj.sigF;
            % Mit bekannter GWL
            else 
                obj.kfs = obj.bestimmekfs(obj.sf,obj.ef,obj.b,obj.c,...
                                      obj.tf,obj.gf,obj.bg,obj.cg,...
                                      obj.sigF,obj.nu,obj.E);
            end
        end
        
        % -----------------------------------------------------------------
        % änderungen für variables k
        if strcmp(obj.kfsopt,'var')
			obj.pwlopt = 'gwl';
            obj.Name = 'PFS_var';
        end
        
        % -----------------------------------------------------------------
        % Bestimme PWL
        if strcmp(obj.pwlopt,'dwl')
            obj.PN = obj.dwl2pwl(obj.E,obj.nu,obj.sf,obj.ef,obj.b,obj.c,obj.ND,...
                                 obj.kfs,obj.nst,obj.miner);
        elseif strcmp(obj.pwlopt,'gwl')
            obj.PN = obj.gwl2pwl(obj.G,obj.tf,obj.gf,obj.bg,obj.cg,obj.ND,...
                                 obj.nst,obj.miner);
        end
        obj.PD = 10^obj.PN(2,end);
        
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
            % ... konstantes k
            if strcmp(obj.kfsopt,'konst') 
                P = gam_a * ( 1 + obj.kfs * snmax );
            % ... variables k (Newtonverfahren)  
            elseif strcmp(obj.kfsopt,'var')
                % ... init mit Startwert
                P = gam_a * ( 1 + obj.kfs * snmax );
                D = obj.damage_akk(P,0);
                if D > 0
                    N = 1/D;
                    %                 N = 1000;
                    f = obj.zielfun(N,gam_a,snmax,obj.gf,obj.tf,obj.bg,obj.cg,...
                        obj.sf,obj.ef,obj.b,obj.c,...
                        obj.sigF,obj.G,obj.nu,obj.nst);
                    % ... Schleife Newtonverfahren
                    %                 fiter = zeros(obj.maxiter,3);   % Debugg speicher
                    iter = 1;
                    while abs(f) > obj.toleranz
                        % ... ableitung (Numerisch)
                        %                     obj.dN = max([floor(N/10000) 1]);
                        %                     fp = obj.zielfun(N+obj.dN,gam_a,snmax,obj.gf,obj.tf,obj.bg,obj.cg,...
                        %                                      obj.sf,obj.ef,obj.b,obj.c,...
                        %                                      obj.sigF,obj.G,obj.nu,obj.nst);
                        %                     fm = obj.zielfun(N-obj.dN,gam_a,snmax,obj.gf,obj.tf,obj.bg,obj.cg,...
                        %                                      obj.sf,obj.ef,obj.b,obj.c,...
                        %                                      obj.sigF,obj.G,obj.nu,obj.nst);
                        %                     dfdN = (fp - fm)/(2*obj.dN);
                        
                        % ... ableitung (analytisch)
                        dfdN = obj.Diffzielfun(N,gam_a,snmax,...
                            obj.gf,obj.tf,obj.bg,obj.cg,...
                            obj.sf,obj.ef,obj.b,obj.c,...
                            obj.sigF,obj.G,obj.nu,obj.nst);
                        % ... Debugg Speicher
                        %                     fiter(iter,:) = [N,f,dfdN];
                        % ... neue Versagensschwingspielzahl
                        N = N - f/dfdN;
                        if N < 1
                            N = 2;
                        elseif N > obj.dauerfest
                            N = obj.dauerfest;
                            break
                        end
                        
                        % ... Zielfunktion
                        f = obj.zielfun(N,gam_a,snmax,obj.gf,obj.tf,obj.bg,obj.cg,...
                            obj.sf,obj.ef,obj.b,obj.c,...
                            obj.sigF,obj.G,obj.nu,obj.nst);
                        
                        % ... keine Endlosschleife
                        iter = iter + 1;
                        if iter > obj.maxiter + 1
                            msg = ['keine Konvergenz für Fatemi Socie Paramter ' ...
                                'iter=',num2str(iter),' ', ...
                                'N=',num2str(N),' ', ...
                                'zielfun=',num2str(f),' '];
%                             warning(msg)
                            break;
                        end
                    end % Ende Newtonverfahren
                    % ... Debugg Speicher
%                     fiter(iter,:) = [N,f,NaN];
                else
                    N = obj.dauerfest;
                end
                % ... Schädigungsparameter (eigentlich quatsch)
                fak = 1;% 
%                 fak = obj.nst^(1/obj.bg);
                P = obj.gf*(2*N*fak)^obj.cg + obj.tf/obj.G*(2*N*fak)^obj.bg;
                P = obj.nst * P;
            end
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
        function [DL,SSP] = lebensdauer(obj,P)
            % Funktion rechnet Lebensdauern aus Schädigungsparametern
            %
            % INPUT:
            % P              - Schädigungsparameter
            %                  1.Zeile Schwingspiele prim. Schädigungsvariable
            %                  2.Zeile Schädigungsparameter
            %                  3.Zeile Durchläufe
            %
            % OUTPUT:
            % DL             - Durchläufe
            % SSP            - Schwingspiele
            %______________________________________________________________
            
            %--------------------------------------------------------------
            % Schädigungsrechnunge
            nP = size(P,2);            % Anzahl Werte in P             
            Dsum = 0;                  % Schadenssumme
            Dlast = 0;                 % Schädigung letzter Durchlauf
            ndl = ceil(max(P(3,:)));   % maximale Anzahl an Durchläufen
            idam = 0;                  % Zeiger auf den Wert an dem Dsum = 1
            ilast = 1;                 % Zeiger auf ersten Wert des letzten Durchlaufs
            for i = 1: nP
                % ... Schädigung aus WL
                [Dakt,Dsum] = obj.damage_akk(P(2,i),Dsum);
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
            
            % Anzahl Punkte in WL
            npt = size(obj.PN,2);
            
            % Sicherstellen P>0
            if P > 0
                % log von P
                logP = log10(P);
                % Dauerfest ?
                if logP > obj.PN(2,npt)
                    % .. suche angrenzente Punkte
                    i = npt;
                    while i>1 && logP>obj.PN(2,i)
                        i=i-1;
                    end
                    % ... Steigung bestimmen N(P)
                    k = ( obj.PN(1,i+1) - obj.PN(1,i) ) / ( obj.PN(2,i+1) - obj.PN(2,i) );
                    
                    % ... Achsabschnitt
                    N0 = obj.PN(1,i) - k * obj.PN(2,i);
                    
                    % ... Linear interpolieren
                    logN = k * logP + N0;
                    
                    % ... Schädigungsbeitrag
                    Dakt = 10^(-logN);
                    
                    % ... Schädigungssumme
                    Dsum = Dsum + Dakt;
                end
                
            end
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
            plot(ax,10.^obj.PN(1,:),10.^obj.PN(2,:))
            set(ax,'XScale','log','YScale','log');
            xlabel('N'), ylabel('P')
            grid on
        end  % Ende plot PWL
    end % Ende Methoden

% -------------------------------------------------------------------------
% STATISCHE METHODEN
% ------------------------------------------------------------------------- 
    methods (Static)
        % ... berechne P-WL aus DWL
        function PN = dwl2pwl(E,nu,sf,ef,b,c,ND,kfs,nst,miner)
            % ... faktor zum verschieben der WL !!! Ausgestellt
            fak = 1;
%             fak = nst^(1/b);            
            % ... WL
            npt = 400;                % Anzahl Punkte auf WL
            PN = zeros(2,npt);      % Speicher (1.Zeile = Ssp, 2.Zeile = P)
            for i = 1:npt
                % ... Punkt N
                N = ND^((i-1)/(npt-1));
                % ... Punkt P
                P = ( (1+nu)*sf/E*(2*N*fak)^b + 1.5*ef*(2*N*fak)^c )...
                        * (1 + kfs * 0.5 * sf * (2*N*fak) ^ b);
                % ... Parallel verschieben WL in Lastrichtung
                P = nst * P;
                % ... Speichern log
                PN(1,i) = log10(N);
                PN(2,i) = log10(P);
            end
            
            % ... Anpassen für elementare oder Mod. Miner
            if miner == 0 || miner == 1
                % ... Steigung letzter Punkt
                k = (PN(1,npt) - PN(1,npt-1))/(PN(2,npt) - PN(2,npt-1));
                if miner == 1 % modifiziert
                    % ... Steigung anpassen bei elementarem Miner
                    k = 2*k - 1;
                end
                % ... neuer Dauerfestigkeitswert bei 1/2
                PN(2,npt) = PN(2,npt) + log10(0.5);
                % ... neue Ecklastspeilzahl bei 1/2 Dauerfestigkeit
                PN(1,npt) = PN(1,npt) + k *log10(0.5);
            end
        end % Ende bestimmen P-WL
        
        % ... berechne P-WL aus GWL
        function PN = gwl2pwl(G,tf,gf,b0,c0,ND,nst,miner) 
            % ... faktor zum verschieben der WL !!! Ausgeschaltet
%             fak =  nst^(1/b0);
            fak = 1;
            % ... WL
            npt = 4;                % Anzahl Punkte auf WL
            PN = zeros(2,npt);      % Speicher (1.Zeile = Ssp, 2.Zeile = P)
            for i = 1:npt
                % ... Punkt N
                N = ND^((i-1)/(npt-1));
                % ... Punkt P  aus Gleitungswl
                P = tf/G*(2*N*fak)^b0 + gf*(2*N*fak)^c0;
                % ... Parallel verschieben WL in Lastrichtung
                P = nst * P;
                % ... Speichern log
                PN(1,i) = log10(N);
                PN(2,i) = log10(P);
            end
            
            % ... Anpassen für elementare oder Mod. Miner
            if miner == 0 || miner == 1
                % ... Steigung letzter Punkt
                k = (PN(1,npt) - PN(1,npt-1))/(PN(2,npt) - PN(2,npt-1));
                if miner == 1 % modifiziert
                    % ... Steigung anpassen bei elementarem Miner
                    k = 2*k - 1;
                end
                % ... neuer Dauerfestigkeitswert bei 1/2
                PN(2,npt) = PN(2,npt) + log10(0.5);
                % ... neue Ecklastspeilzahl bei 1/2 Dauerfestigkeit
                PN(1,npt) = PN(1,npt) + k *log10(0.5);
            end
        end % Ende bestimmen P-WL
        
        % ... bestimme kfs parameter
        function [kds,k] = bestimmekfs(sf,ef,b,c,tf,gf,b0,c0,sigF,nu,E)
            % Berechnet k für Fatemi Socie Parameter
            %
            % INPUT:
            % sf,ef,b,c     - Parameter DWL
            % tf,gf,b0,c0   - Parameter GWL
            % sigF          - Fließspannung
            % nu            - Querdehnzahl (elastisch)
            %
            % OUTPUT:
            % kds           - k/sigF
            % k             - k Parameter Fatemi Socie Parameter
            %
            % ANMERKUNG:
            % Vorgehen aus D. McClaflin, A. Fatemi 2003 - Torsional deformation
            %              and fatigue of hardened steel including mean stress
            %              and stress gradient effects
            %
            % Bestimme k für N = 10^2...10^7 und nehme Mittelwert
            %
            %______________________________________________________________
            
            % ... berechne Schubmodul
            G = E/(2*(1+nu));
            
            % ... Bestimme Punkte und k über Mittelwert (log Verteilt)
%             Nmin = 3;
%             Nmax = 7;
%             ndp = Nmax-Nmin+1;
%             N = logspace(Nmin,Nmax,ndp);
%             fak = ((tf/G*(2*N).^b0 + gf*(2*N).^c0) ./ ...
%                 ( (1+nu)*sf/E*(2*N).^b + 1.5*ef*(2*N).^c ) - 1);
%             k = fak * ...
%                 2*sigF./(sf * (2*N).^b);
%             k = 1/length(k) * sum(k);
            
            
            % ... bestimme Integralen Mittelwert über Simpson Integration
            N1 = 1e3; N2 = 1e7;
            nip = 1001;         % Anzahl Inkremente + 1
            h = (N2-N1)/(nip-1);
            N = N1:h:N2;
            ksimpson = ((tf/G*(2*N).^b0 + gf*(2*N).^c0) ./ ...
                ((1+nu)*sf/E*(2*N).^b + 1.5*ef*(2*N).^c ) - 1)* ...
                2*sigF./(sf * (2*N).^b);
            F = ksimpson(1) + ksimpson(nip);
            for i = 2:nip-1
                if mod(i,2) == 0
                    F = F + 2*ksimpson(i);
                else
                    F = F + 4*ksimpson(i);
                end
            end
            k = h/3 * F/(N2-N1);
            
            % ... k/sigF
            kds = k/sigF;
            
        end % Ende bestimme kfs
        
        % ... Zielfunktion
        function error = zielfun(N,gam_a,snmax,...
                gf,tf,b0,c0,...
                sf,ef,b,c,...
                sigF,G,nu,nst)
            % !! k als Werkstoffkennwert ohne Stuetzwirkung
            % ... Stützwirkung
%             fak = nst^(1/b0);
            fak = 1;% 
            % ... GleitungsWL
            gwl = gf*(2*N*fak)^c0 + tf/G*(2*N*fak)^b0;
%             gwl = nst * gwl;
            % ... angepasse DehnungsWL
            E = G * 2 * (1+nu);
            fak = 1; % 
%             fak = nst^(1/b);
            ewl = (1+nu)*sf/E*(2*N*fak)^b+1.5*ef*(2*N*fak)^c;
%             ewl = nst * ewl;
            % ... k = (gwl/ewl -1) * sigF/(sf*(2*N)^b)
            k = (gwl/ewl-1)*2*sigF/(sf*(2*N*fak)^b);
            % ... PFS
            PFS = gam_a*(1+k*snmax/sigF);
            % ... zielfunktion
            error = PFS - gwl;
        end % Ende Zielfunktion       
        
        % ... Ableitung der Zielfunktion
        function Diff = Diffzielfun(N,gam_a,snmax,...
                gf,tf,b0,c0,...
                sf,ef,b,c,...
                sigF,G,nu,nst)
            % ... Stützwirkung
            fakg = nst^(1/b0);
%             fakg = 1;% 
            % ... GleitungsWL
            gwl = gf*(2*N*fakg)^c0 + tf/G*(2*N*fakg)^b0;
%             gwl = nst * gwl;
            % ... angepasse DehnungsWL
            E = G * 2 * (1+nu);
            fake = nst^(1/b);
%             fake = 1;
            ewl = (1+nu)*sf/E*(2*N*fake)^b+1.5*ef*(2*N*fake)^c;
%             ewl = nst * ewl;
            % ... Ableitung GleitungsWL
            Diffgwl = (2*fakg)^b0*b0*tf/G*(N)^(b0-1) + (2*fake)^c0*c0*gf*(N)^(c0-1);
%             Diffgwl = nst * Diffgwl;
            % ... Ableitung DehnungsWL
            Diffewl = (1+nu)*sf/E*(2*fake)^b*b*(N)^(b-1) + 3/2*ef*(2*fake)^c*c*(N)^(c-1);
%             Diffewl = nst * Diffewl;
            % ... Ableitung PFS
            DiffPFS = gam_a*snmax/sigF*2*sigF/sf*(...
                      -(2*fake)^(-b)*b*(N)^(-b-1) * (gwl/ewl-1) ...
                      + (2*fake*N)^(-b) * (Diffgwl*ewl-gwl*Diffewl)/ewl^2);
            % ... Ableitung
            Diff = DiffPFS - Diffgwl;
        end % Ende Zielfunktion  
        
        % ... Interpolation zwischen Gestaltänderungsergie- und
        % Normaldehnungshypothese
        function q = skaliereGEHundNH(fwt)
            % q(x) = a - b * f(x)
            % q(1/sqrt(3)) = 0    -> GEH
            % q(1)         = 1    -> NH
            % f(1/sqrt(3)) = f13
            % f(1)         = f1
            % a =  f13/(f13-f1)
            % b = 1/(f13-f1)
            % Hier Ansätze  f(x) = (1/x)^d (d=1 entspricht RiLiNiLi)
            % -------------------------------------------------------------
            d = 1;
            f1 = 1;
            f13 = sqrt(3)^d;
            q = (f13 - 1./fwt.^d)./(f13-f1);
        end
        
    end % Ende statische Methoden

end % Ende Klassendefinition PRAM
