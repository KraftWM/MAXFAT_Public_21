classdef PZ2_stuetz < handle
% -------------------------------------------------------------------------
% Klassendefinition: Schädigungsparameter PZ mit
% Wöhlerliniendefinition aus Stützstellen
% Zusammenfassung aller Definitionen/ Optionen & Funktionen für
% Schädigungsparameter PZ, nach Hertel
%
%
%__________________________________________________________________________
% KONSTRUKTOR:
%    P = PZ(E,nu,sigF,Weff_WS_stuetz,d,fak,varargin);
%           E             - E Modul
%           nu            - Querdehnzahl
%           sigF          - Fließspannung
%          Weff_WS_stuetz - Stützstelle Weff-Wöhlerlinie
%             d           - Neigung
%            fak          - Faktor Pz_bauteil = 1/fak * Pz_WS_stuetz
%                          Zusammenfassung von Sicherheits & Stützfaktoren
%           varargin      - variabler Input
%                       'VarName1',var1,'VarName2',var2,...)
% mit varargin können alle anderen zusätzlichen public Eigenschaften
% gesetzt werden
%__________________________________________________________________________
%
% EIGENSCHAFTEN:
% Public:
%   E            - E Modul
%   nu           - Querdehnzahl
% Weff_WS_stuetz - Stützstelle Weff-Wöhlerlinie
% Weff_WSD_stuetz- Stützstelle Weff-Wöhlerlinie Dauerfestigkeit
%    d           - Neigung
%   fak          - Faktor Pz_bauteil = 1/fak * Pz_WS_stuetz
%                          Zusammenfassung von Sicherheits & Stützfaktoren
% ce             - Endrisslänge (auf der Oberfläche, Punkt B)
% Gsig,Gtau      - Bezogener Spannungsgradient
%    M_sig       - Mittelspannungsempfindlichkeit (nach FKM RiLiNiLi aus
%                    Zugfestigkeit) % !!!!! vlt noch zu modifizieren  
% ators          - Einfluss der Torsionsspannung auf Rissschließen
%  mu = 0.5;     - Reibkontakt der Rissufer, Einflus sn
% takt = 50;     - Reibkontakt der Rissufer, zu überwindende Schubspannung
%  sigF          - Fließspannung (aus Rm und Rp02)
% fwt            - Schubwechselfestigkeit zu Wechselfestigkeit
%                    default fwt = 1/sqrt(3) -> Mises
% tauF           - Fließspannung Torsion
% FA = 1000      - Abklingkonstante Rissöffnung      
%  N_stuetz        - Stützstelle Lebensdauer
%                    default N_stuetz = 1
%  ND_stuetz       - Stützstelle Lebensdauer Dauerfestigkeit
%                    default N_stuetz = 1e6
%
% Protected:
%  InterpMethod     - Methode für Interpolation der Gewichtsfunktion (siehe
%                     Matlab Funktion griddedInterpolant, default = 'linear')
%  YIA              - griddedInterpolant für Geometriefunktionen
%  YIIA             - griddedInterpolant für Geometriefunktionen
%  YIIIA            - griddedInterpolant für Geometriefunktionen
%  YIB              - griddedInterpolant für Geometriefunktionen
%  YIIB             - griddedInterpolant für Geometriefunktionen
%  YIIIB            - griddedInterpolant für Geometriefunktionen
%   PN              - Logaritmierte Werte der P-WL
%    a0             - Anfangsrisstiefe
%    c0             - Anfangsrissläng
%    ae             - Endrisslänge tiefe
%    ce             - Endrisslänge oberfläche
%   a0zc0           - Verhältniss a0 zu c0
%   Q               - Stützpunkt PWL
%   PZD0            - Dauerfestigkeit PZ
%   ND0             - Abknickpunkt dauerfestigkeit
%   mJ              - Steigung der PZ-BasisWöhlerlinie
%   CJ              - Rissfortschrittskonstante
%   ZeffthLR        - Langrissschwllenwert 
%   lstern          - mikrostrukturelle Hilfsgröße
%   Y1A0            - Geometriefunktion Mode I Punkt A,   Startwert
%   Y2A0            - Geometriefunktion Mode II Punkt A,  Startwert
%   Y3A0            - Geometriefunktion Mode III Punkt A, Startwert
% exmax             - Bisher aufgetretene Maximal(normal)dehnung 
% exmin             - Bisher aufgetretene Minimal(normal)dehnung 
% exop_alt          - Rissöffnungs(normal)dehnung aus Vorgeschichte
% exop_ein          - Rissöffnungs(normal)dehnung für einstufige Belastung 
%                    (aus letztem Ssp)
% speicherrifo      - Speicher für Rissfortschritt
%                     ['DL','ZYK','ai[µm]','ci[µm]','ai/ci']
%
% Name              - Name des Parameters (nur zum generieren von Output)
%__________________________________________________________________________
% FUNKTIONEN:
% Normal:
%      dwl2weffwl     -> berechne Weff - Wöhlerlinie aus Dehnungswöhlerlinie
%      rainflow       -> rainflowzählung mit hcm
%    kurzriss_mode1   -> berechne PZ bei Mode I
%    kurzriss_mode23  -> berechne PZ bei Mode II & III
%    lebensdauer      -> Lebensdauerrechnung
%    inkrifo          -> inkrementeller Rissfortschritt
%    plotPWL          -> Plotet PWL
% Static:
%      closesig       -> Rissschließspannung aus Rissschließdehnung für
%                        Ramberg Osgood Gleichung 
%      sigclose       -> Rissschließspannung aus Rissschließdehnung für
%                        beliebige Hystereseäste
%   Y_OFR_TYPB_ZD     -> Geometriefunktionen Oberflächenriss Typ B Riss ZD
%   Y_OFR_TYPA_ZD     -> Geometriefunktionen Oberflächenriss Typ A Riss ZD
%   Y_OFR_TYPA_SH     -> Geometriefunktionen Oberflächenriss Typ B Riss SH
%  mises_vorzeichen   -> Berechnet die Vorzeichenbehaftete von Mises 
%                        Vergleichsspannungen
%  newman_max         -> bestimme Rissöffnungsspannung
%  epsopen            -> bestimme Rissöffnungsdehnung aus Rissöffnungsspannung
% verzerrungsenergie  -> berechnen der Verzerrungsenergie der Hysteresehalbäste 
%                        Integration mit Trapezregel
% DamAkk              -> Schadensakkumulation
% simpson             -> simpson integration
% skaliereGEHundNH    -> Interpolation zwischen GEH und NH
%__________________________________________________________________________
% EXTRERNE FUNKTIONEN:
% akima2d       -> 2D interpolation mit splines (schneller als matlab
%                  version)
% polyfit       -> Matlab funktion um Daten mit polynom zu fitten
%__________________________________________________________________________
% -------------------------------------------------------------------------
% AUTHOR: JAN KRAFT 
% STAND: September 2021
% -------------------------------------------------------------------------  

% -------------------------------------------------------------------------
% EIGENSCHAFTEN
% -------------------------------------------------------------------------
% PUBLIC:
   properties % Notwendige Größen für Kerbnäherung
       E {mustBeNumeric} = NaN;                                            % E Modul
       nu {mustBeNumeric} = NaN;                                           % Querdehnzahl
       Pz_WS_stuetz {mustBeNumeric} = NaN;
       Pz_WSD_stuetz {mustBeNumeric} = NaN; 
       d {mustBeNumeric} = NaN;
       fak   {mustBeNumeric} = 1;    
       N_stuetz {mustBeNumeric} = 1;   
       ND_stuetz = 1e6;
       ce {mustBeNumeric} = 0.25;                                          % Endrisslänge (auf der Oberfläche, Punkt B)
       Gsig {mustBeNumeric} = 0;                                           % 
       Gtau {mustBeNumeric} = 0;                                           % Bezogener Spannungsgradient
       M_sig {mustBeNumeric} = -1;                                         % Mittelspannungsempfindlichkeit (nach FKM RiLiNiLi aus
                                                                           % Zugfestigkeit) % !!!!! vlt noch zu modifizieren
       ators {mustBeNumeric} = 0;                                          % Einfluss torsionsspannung auf rissschließen                           
       sigF {mustBeNumeric} = NaN;                                         % Fließspannung (aus Rm und Rp02)
       tauF {mustBeNumeric} = NaN;                                         % Fließspannung Schub
       mu {mustBeNumeric} = 0.5;                                           % Reibkontakt der Rissufer, Einflus sn
       takt {mustBeNumeric} = 50;                                          % Reibkontakt der Rissufer, zu überwindende Schubspannung
       FA {mustBeNumeric} = 1000;                                          % Abklingkonstante Rissöffnung                                                             
       exmax {mustBeNumeric} = 0;                                          % Bisher aufgetretene Maximaldehnung
       exmin {mustBeNumeric} = 0;                                          % Bisher aufgetretene Minimaldehnung
       exop_alt {mustBeNumeric} = 0;                                       % Rissöffnungsdehnung aus Vorgeschichte
       exop_ein {mustBeNumeric} = 0;                                       % Rissöffnungsdehnung für einstufige Belastung (aus letztem Ssp)

       fwt {mustBeNumeric} = 1/sqrt(3);                                    % Schubwechselfestigkeit zu Wechselfestigkeit
       
   end % Ende EIGENSCHAFTEN Public 
   
% PROTECTED
   properties (SetAccess = protected) % Größen werden automatisch gesetzt
       InterpMethod = 'linear';                                            % Interpolationsvorschrift Gewichtsfunktionen
       YIA                                                                 % griddedInterpolant für Geometriefunktionen
       YIIA                                                                % griddedInterpolant für Geometriefunktionen
       YIIIA                                                               % griddedInterpolant für Geometriefunktionen
       YIB                                                                 % griddedInterpolant für Geometriefunktionen
       YIIB                                                                % griddedInterpolant für Geometriefunktionen
       YIIIB                                                               % griddedInterpolant für Geometriefunktionen
       PN {mustBeNumeric} = NaN;                                           % P-WL
       dauerfest = 1e20;                                                   % definition für dauerfest
       a0 {mustBeNumeric} = NaN;                                           % Anfangsrisstiefe
       c0 {mustBeNumeric} = NaN;                                           % Anfangsrissläng
       ae {mustBeNumeric} = NaN;                                           % Endrisslänge tiefe
       a0zc0 {mustBeNumeric} = NaN;                                        % Verhältniss a0 zu c0
       Q {mustBeNumeric} = NaN;                                            % Stützpunkt Wöhlerlinie
       PZD0 {mustBeNumeric} = NaN;                                         % Dauerfestigkeit PZ
       ND0 {mustBeNumeric} = NaN;                                          % Abknickpunkt Dauerfestigkeit  
       mJ {mustBeNumeric} = NaN;                                           % Steigung der PZ-BasisWöhlerlinie
       CJ {mustBeNumeric} = NaN;                                           % Rissfortschrittskonstante
       ZeffthLR {mustBeNumeric} = NaN;                                     % Langrissschwllenwert 
       lstern {mustBeNumeric} = NaN;                                       % mikrostrukturelle Hilfsgröße
       Y1A0 {mustBeNumeric} = NaN;                                         % Geometriefunktion Mode I Punkt A
       Y2A0 {mustBeNumeric} = NaN;                                         % Geometriefunktion Mode II Punkt A
       Y3A0 {mustBeNumeric} = NaN;                                         % Geometriefunktion Mode III Punkt A
       speicherrifo {mustBeNumeric} = NaN;                                 % Speichert Risstiefe (= a am Punkt A, nur für auswertung)
       Name = 'PZ_st';                                                     % Name des Parameters
       Pz_Bau_stuetz {mustBeNumeric} = NaN;                                % Stuetzstelle Bauteil
       Pz_BauD_stuetz {mustBeNumeric} = NaN;                               % Stuetzstelle Bauteil
   end % Ende EIGENSCHAFTEN protected 
   

% -------------------------------------------------------------------------
% KONSTRUKTOR
% ------------------------------------------------------------------------- 
methods
    function obj = PZ2_stuetz(E,nu,sigF,Pz_WS_stuetz,d,fak,varargin)
        %------------------------------------------------------------------
        % Konstruktor der Klasse
        % INPUT:
        %      E,nu,             - Elastizität
        %      sigF              - Fließspannung
        %      Pz_WS_stuetz      - Stützstellen Pram-Wöhlerlinie bei definierter
        %                          Lebensdauer für Werkstoff
        %                          Pz_stuetz bei N = 1
        %         d              - Neigungen PWL
        %                          d1 für N < 1000
        %           varargin     - variabler Input
        %                       'VarName1',var1,'VarName2',var2,...)
        % OUTPUT:
        %    obj - objekt der Klasse
        %------------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % Setze Felder
        obj.E = E;
        obj.nu = nu;
        obj.sigF = sigF;
        obj.Pz_WS_stuetz = Pz_WS_stuetz;
        obj.d = d;
        obj.fak = fak;
        
        % -----------------------------------------------------------------
        % Variabler Input
        nvarargin = length(varargin);
        if nvarargin > 0
            for i = 2 : 2 : nvarargin
                varname = varargin{i-1};
                varvalue = varargin{i};
                switch varname
                    case 'ce'
                        obj.ce = varvalue;
                    case 'M_sig'
                        obj.M_sig = varvalue;
                    case 'Gsig'
                        obj.Gsig = varvalue;
                    case 'Gtau'
                        obj.Gtau = varvalue;
                    case 'FA'
                        obj.FA = varvalue;
                    case 'mu'
                        obj.mu = varvalue;
                    case 'takt'
                        obj.takt = varvalue;
                    case 'ators'
                        obj.ators = varvalue;
                    case 'fwt'
                        obj.fwt = varvalue;
                    case 'tauF'
                        obj.tauF = tauF;
                    case 'Pz_WSD_stuetz'
                        obj.Pz_WSD_stuetz = varvalue;
                    case 'N_stuetz'
                        obj.N_stuetz = varvalue;
                    case 'ND_stuetz'
                        obj.ND_stuetz = varvalue;
                    otherwise
                        msg = ['Variable ', varname, ' nicht erkannt'];
                        warning(msg)                        
                end % Ende switch variablenname
            end % Ende Schleife über variablen input
        end
        
        % -----------------------------------------------------------------
        % Interpolationsobjekte für Geometriefunktionen
        [~,obj.YIIB,obj.YIIIA,~] = obj.Y_OFR_TYPA_SH();
        [obj.YIA,obj.YIB] = obj.Y_OFR_TYPA_ZD();
        [~,~,obj.YIIA,~,~,obj.YIIIB] = obj.Y_OFR_TYPB_ZD();
        
            
        % -----------------------------------------------------------------
        % Grenzen fwt
        if obj.fwt < 1/sqrt(3)   % Untere Grenze von Mises
            obj.fwt = 1/sqrt(3);
        elseif obj.fwt > 1       % Obere Grenze Normalenhypothese
            obj.fwt = 1;
        end

        % -----------------------------------------------------------------
        % Fließspannung Torsion
        if isnan(obj.tauF)
            q = obj.skaliereGEHundNH(obj.fwt);
            tfGEH = obj.sigF/sqrt(3);
            tfNH = obj.sigF;
            obj.tauF = q * tfNH + (1-q) * tfGEH;
        end
        
        % -----------------------------------------------------------------
        % bestimme P-Wöhlerlinie
        mJ = -1/d;        
        Q = obj.Pz_WS_stuetz^mJ;
%         ndp = 21;  % Anzahl Datenpunkte auf PWL
%         [mJ,Q,~,Weff] = obj.dwl2weffwl(ndp);
        
        % ------------------------------------------------------------------
        % a/c Verhältnis bestimmen sodass a/c = konst, d(a/c)/dn = 0
        % ... Init a/c
        azc = 1;
        % ... Geometriefaktoren
%         [Y1A,Y1B] = obj.Y_OFR_TYPA_ZD(azc,0);
        Y1A = obj.YIA(azc,0);
        Y1B = obj.YIB(azc,0);
        % ... Iteration a/c
        f = abs((Y1A/Y1B)^(2*mJ) - azc);
        while f > 1e-4
            azc = 0.5 * ( azc + (Y1A/Y1B)^(2*mJ)); % Neuer Wert = mittelwert aus alt+a/c aus Rissfortschritt
%             [Y1A,Y1B] = obj.Y_OFR_TYPA_ZD(azc,0);
            Y1A = obj.YIA(azc,0);
            Y1B = obj.YIB(azc,0);
            f = abs((Y1A/Y1B)^(2*mJ) - azc);
        end
        % Hier ist Bisektion vielleicht stabiler
        
        % -----------------------------------------------------------------
        % Q bestimmen mit Geometriefunktion (Lebensdauer bei PZ = 1)
%         Q = Q * (2*pi*Y1A^2)^mJ;
        
        % -----------------------------------------------------------------
        % Kurzrissinitiierungsschwellwert
        % ... Nicht Selbst definiert bei N_Stütz festlegen
        if isnan(obj.Pz_WSD_stuetz)
            obj.Pz_WSD_stuetz = obj.Pz_WS_stuetz *(obj.ND_stuetz/obj.N_stuetz)^(obj.d);
        end
        % ... PZ Dauerfestigkeit
%         PZD0 = 2*pi*Y1A^2*obj.Pz_WSD_stuetz;
        PZD0 = obj.Pz_WSD_stuetz;
        % ... Neuer Abknickpunkt
        ND0 = Q * PZD0^(-mJ);
        
        % -----------------------------------------------------------------
        % CJ und Langrissschwellwert nach Vormwald abschätzen 
        % ... Rissfortschritskonstante
        CJ = 10^(-5)*(5*10^5/E)^mJ;
        % ... Langrissschwellwert
        ZeffthLR = E/(5*10^6);

        % -----------------------------------------------------------------
        % Rissfortschritt rückwärts zur Ermittlung von a0 und c0 bei PZ=1
        % ... Endrisslänge Punkt A
        ae = obj.ce * azc;
        % ... Anfangsrisslänge
        a0 = ( ae^(1-mJ) - (1-mJ)*CJ*Q )^(1/(1-mJ));
%         c0 = a0/azc;
        
        % -----------------------------------------------------------------
        % mikrostrukturelle Hilfsgröße bzw. Zusatzrisslänge l* bzw. Größe 
        % der plastischen Zone am gestoppten Riss
        lstern = ZeffthLR/PZD0 - a0;
        
        % -----------------------------------------------------------------
        % Geometriefunktionen 
        % (kurzer Riss a0/c0, ungekerbt a0/rho=0, tiefster Punkt A)
        % ... Mode I
        Y1A0 = Y1A;
        % ... Mode II
%         [~,~,Y2A0] = obj.Y_OFR_TYPB_ZD(azc,0);
        Y2A0 = obj.YIIA(azc,0);
        Y2A0 = 2*Y2A0; % Wegen Hinweis Hertel zur FE Lösung
        % ... Mode III
%         [~,~,Y3A0] = obj.Y_OFR_TYPA_SH(azc,0);
        Y3A0 = obj.YIIIA(azc,0);
        
        % -----------------------------------------------------------------
        % Stützwirkung
        % ... Verschiebung der PZ-WL 
%         Q = Q * obj.nst^(-1/obj.bf);
%         Q = Q * (1/obj.fak)^(2*mJ);
        Q = Q * (1/obj.fak)^(mJ);
        % ... neues a0
        a0 = ( ae^(1-mJ) - (1-mJ)*CJ*Q )^(1/(1-mJ));
        % ... neues c0
        c0 = a0/azc;
        % ... neue PZ Dauerfestigkeit
        PZD0 = 1/obj.fak * obj.Pz_WSD_stuetz; % ZeffthLR/(a0+lstern);
        % ... neuer Abknickpunkt
        ND0 = Q * PZD0^(-mJ);
        
        % -------------------------------------------------------------------------
        % Zusammenfassung aller restlichen Kurzrissvariablen
        % ... Rissgeometrie
        obj.ae = ae;                                                       % Endrisstiefe [mm]
        obj.a0 = a0;                                                       % Anfangsrisstiefe [mm]
        obj.c0 = c0;                                                       % (halbe) Anfangsrisslänge [mm]
        obj.a0zc0 = azc;                                                   % Risslängenverhältniss [mm]
        obj.Y1A0 = Y1A0;                                                   % Geometriefaktor Mode I
        obj.Y2A0 = Y2A0;                                                   % Geometriefaktor Mode II
        obj.Y3A0 = Y3A0;                                                   % Geometriefaktor Mode III
        % ... Rissfortschrittsparameter
        obj.mJ = mJ;                                                       % Exponent Rissfortsgleichung [-]
        obj.CJ = CJ;                                                       % [mm/cyc*(N/mm)^mJ]
        obj.Q = Q;
        obj.PZD0 = PZD0;
        obj.ND0 = ND0;
        obj.lstern = lstern;                                               % [mm]
        obj.ZeffthLR = ZeffthLR;       
    end % Ende Konstruktor
end % Ende Konstruktor

% -------------------------------------------------------------------------
% METHODEN
% ------------------------------------------------------------------------- 
methods
%     % ... bestimme Weff Wöhlerlinie
%     function [mJ,Q,N,Weff] = dwl2weffwl(obj,ndp)
%         % -----------------------------------------------------------------
%         % Funktion bestimmt Weff Wöhlerlinie aus DWL(glatte Probe R=-1)
%         % INPUT:
%         % ndp           - Anzahl Datenpunkte auf WL
%         % OUTPUT:
%         % mJ            - Rissfortschrittsparameter (Steigung WL)
%         % Q             - Abknickpunkt PZ-WL (bei PZ=1, hier erstmal Weff=1)
%         % N             - Stützpunkte Lebensdauern
%         % Weff          - Stützpunkte Verzerrungsenergie
%         % -----------------------------------------------------------------
%         % ... Ramberg Osgood Parameter aus Konsostenzbedingung
%         ns = obj.bf/obj.cf;
%         Ks = obj.sf/obj.ef^ns;
%         
% %         ns = obj.nstrich;
% %         Ks = obj.Kstrich;
%         
%         % ... Speicher
%         N = zeros(1,ndp);   % Speicher Lebensdauer Wöhlerlinie
%         Weff = zeros(1,ndp);% Speicher effektive Verzerrungserngiedichte WL
%         
%         % ... Datenpunkte erzeugen
%         for i = 1 : ndp
%             % ... Lebensdauern von 100 bis ND
%             N(i) = 100 * (obj.ND/100)^((i-1)/(ndp-1));
% 
%             % ... Spannungs- und Dehnungsamplitude aus MBC Gleichung (ohne nst)
%             sa = obj.sf*(2*N(i))^obj.bf;
%             ea = sa/obj.E + obj.ef*(2*N(i))^obj.cf; 
% 
%             % ... Spannungs- und Dehnungsamplitude mit stat. Stützwirkung
%         % 	sa =        nst           * sf*(2*N(i))^b;   % Verschieben sig-WL mit
%         % 	ea = sa/E + nst^(1/ns) * ef * (2*N(i))^c;   % Beachtung der Kompatibilität!!!
% 
%             % ... Rissöffnungsspannung
%             sop = sa*( 0.535*cos(pi/2*sa/obj.sigF) - 0.344* sa/obj.sigF);
%             sop = max( [-sa min( [sa, sop] )] );
% 
%             % ... Rissschließdehnung
%             ecl = -ea + (sop+sa)/obj.E + 2 * ((sa+sop)/(2*Ks))^(1/ns);
% 
%             % ... Rissschließspannung
%             scl = obj.closesig(sa,ea,ecl,obj.E,Ks,ns);
%             scl = max( [-sa min( [sa scl] )] );
% 
%             % ... effektive Verzerrungsenergiedichte
% %             Weffe(i) = (sa-scl)^2/(2*obj.E);
% %             Weffp(i) = (sa-scl)*((ea-ecl)-(sa-scl)/obj.E)/(1+ns);
%             Weff(i) =  (sa-scl)^2/(2*obj.E) + (sa-scl)*((ea-ecl)-(sa-scl)/obj.E)/(1+ns);
% 
%         end % Ende Datenpunkte WL
% 
% 
%         % -----------------------------------------------------------------
%         % Steigung der Weff WL aus lin. Regression der log. Werte;
%         % ... log Werte
%         N = log10(N);
%         Weff = log10(Weff);
%         % ... lineare Regression
%         p = polyfit(Weff,N,1);
%         % ... Auslesen Werte
%         mJ = -p(1);                      % Steigung
%         Q = 10^p(2);                     % Lebendsdauer bei Weff = 1;
%     end % Ende Bestimmen Weff Wöhlerlinie
%     
    % ... Rainflowzählung
    function P = rainflow(obj,sigepsfile,ndata,phi,psi)
        % HCM Hysteresis Counting Method für Kurzrissmodell
        %
        % Im Residuum werden keine Werte sondern nur Spaltenindices
        %       gespeichert
        %
        % Der vorliegende Code wurde auf Grundlage von
        % RAINFLOW-HCM Ein Hyseresisschleifen-Zählalgorithmus auf
        % Werkstoffmechanischer Grundlage von U.H. Chlormann und T. Seeger aus dem
        % Jahr 1985 implementiert
        %
        %==================================================================
        % INPUT:
        %
        %
        % Aus Kerbsimulation
        % sigepsfile        - Verweis auf Datei mit diskreten
        %                     Lastpunkten für (elast-plast) Spannungen
        %                     und Dehungen im Kerbkoordinantensystem
        %                     Data = [DLZ,S11,S22,S12,E11,E22,E33,G12]
        %                     wird automatisch umgewandelt in
        %                     Data = [S,E,DLZ];
        %                     (double array e R^(13 x numdata))
        %                     DLZ -> Durchlaufzähler
        % ndata             - Anzahl Zeitpunkte in Lastfolge
        % phi,psi           - Winkel kritische Ebene (in Rad)
        %
        % Geschichtsvariablen Kurzrissmodell:
        % exmax             - Bisher aufgetretene Maximaldehnung
        % exmin             - Bisher aufgetretene Minimaldehnung
        % exop_alt          - Rissöffnungsdehnung aus Vorgeschichte
        % exop_ein          - Rissöffnungsdehnung für einstufige Belastung
        %                    (aus letztem Ssp)
        %
        % Geometriefunktionen Kurzriss
        % Y1A0              - Geometriefaktor Riss Mode I am Punkt A
        % Y2A0              - Geometriefaktor Riss Mode II am Punkt A
        % Y3A0              - Geometriefaktor Riss Mode III am Punkt A
        %
        % Materialparameter:
        % FA                - Abklingfaktor Rissöffnungsdehnung (nach Vormwald)
        % Q                 - Lebensdauer bei PZ=1
        % mJ                - Steigung der PZ-Wöhlerlinie
        % sf                - Fließspannung
        % takt              - Schubspannung für Reibkontakt der Rissufer
        % mu                - Faktor zur Berücksichtigung des Einfluss der
        %                     Normalspannung auf Rissuferkontakt
        % nu                - Querdehnzahl
        %
        % OUTPUT:
        % P                 - Schädingungsparameter
        %                      (1.Zeile) -> Rissöffnungsmode
        %                      (2.Zeile) -> gezählte Schwingspiele
        %                      (3.Zeile) -> PZ Parameter
        %                      (4.Zeile) -> Durchlauf in dem Ssp gezählt wurde
        %==================================================================
        % Autor Jan Kraft
        % Version 4, April 2021
        %==================================================================
        
        %==================================================================
        % Initialisierung
        % ... ein paar Definitionen
        nBuf = 10000;                                                      % Buffergröße
        nBufMax = 100e6;                                                   % Buffergröße die nicht überschritten werden soll
        nSp = 13;                                                          % Anzahl der Spalten (1DLZ 3Spannungen 4 Dehnungen)
        datacounter = 0;                                                   % Zählt wieviele Daten gelesen wurden
        % ... Öfnne Lastdatei & Lese Header
        fidLoad = fopen(sigepsfile,'r');                                   % FileID
        % ... Vorbereiten Rainflow
        nRes = min(ndata,5e5);                                             % Beschränke Größe des Residuums
        nDam = min(ndata,5e5);                                             % Beschränke Größe der Outputvariable
        IR1 = 1;                                                           % Zeiger auf den letzten nicht schließfähigen Wert (Mode I)
        IR2 = 1;                                                           % (Mode II)
        IR3 = 1;                                                           % (Mode III)
        IZ1 = 1;                                                           % Zeiger auf den letzten noch schließfähigen Wert (Mode I)
        IZ2 = 1;                                                           % (Mode II)
        IZ3 = 1;                                                           % (Mode III)
        counter1 = 0;                                                      % Schleifenzähler (Mode I)
        counter2 = 0;                                                      % (Mode II)
        counter3 = 0;                                                      % (Mode III)
        RES1 = zeros(1,nRes); RES1(1) = 1;                                 % Residuum (Mode I)
        RES2 = zeros(1,nRes); RES2(1) = 1;                                 % (Mode II)
        RES3 = zeros(1,nRes); RES3(1) = 1;                                 % (Mode III)
        Buf1 = zeros(nSp,nBuf);                                            % Speicher (Mode I)
        Buf2 = zeros(nSp,nBuf);                                            % (Mode II)
        Buf3 = zeros(nSp,nBuf);                                            % (Mode III)
        nextDat1 = 1;                                                      % Zeiger auf nächste freie Stelle im Buffer (Mode I)
        nextDat2 = 1;                                                      % (Mode II)
        nextDat3 = 1;                                                      % (Mode III)
        empty1 = nBuf;                                                     % Anzahl freie Stellen im Buffer (Mode I)
        empty2 = nBuf;                                                     % (Mode II)
        empty3 = nBuf;                                                     % (Mode III)
        Nf1 = 2:nBuf+1; Nf1(nBuf) = 1;                                     % Zeiger auf Nachfolger im Buffer (Mode I)
        Nf2 = 2:nBuf+1; Nf2(nBuf) = 1;                                     % (Mode II)
        Nf3 = 2:nBuf+1; Nf3(nBuf) = 1;                                     % (Mode III)
        Vg1 = 0:nBuf-1; Vg1(1) = nBuf;                                     % Zeiger auf Vorgänder im Buffer (Mode I)
        Vg2 = 0:nBuf-1; Vg2(1) = nBuf;                                     % (Mode II)
        Vg3 = 0:nBuf-1; Vg3(1) = nBuf;                                     % (Mode III)
        % ... Init Zeiger auf freie stelle in Schädigung
        P = zeros(4,nDam);                                                 % Speicher für Schädigungsparameter
        p = 0;                                                             % zeiger auf nächste frei Spalte in P
        % ... Rotationsmatricen von Spannungen und Dehnungen in
        % kritische Ebene zu Drehen
        PS = transformstress(phi,psi);                                     % Dreht Spannungen
        PE = transformstrain(phi,psi);                                     % Dreht Dehnungen
        
        %==============================================================
        % Rainflow Zählung
        % Äußere Schleife über alle Werte in sigepsfile
        status = 1;                                                    % Abbruchbedingung
        while status && empty1 > 1 && empty2 > 1 && empty3 > 1
            idx1 = NaN(1,nBuf);                                            % init Indices in Buffer in die beschrieben werden (Mode I)
            idx2 = NaN(1,nBuf);                                            % (Mode II)
            idx3 = NaN(1,nBuf);                                            % (Mode III)
            ndatared = 0;                                                  % Anzahl (in diesem Durchlauf) gelesener Daten
            % 1. innere Schleife zum füllen des Buffers
            while status && ...
                  Nf1(nextDat1) ~= RES1(IR1) && ...
                  Nf2(nextDat2) ~= RES2(IR2) && ...
                  Nf3(nextDat3) ~= RES3(IR3)
                % merke Index
                idx1(ndatared+1) = nextDat1;
                idx2(ndatared+1) = nextDat2;
                idx3(ndatared+1) = nextDat3;
                % Zeiger auf nächste freie Stelle Setzten
                if datacounter < ndata
                    nextDat1 = Nf1(nextDat1);
                    nextDat2 = Nf2(nextDat2);
                    nextDat3 = Nf3(nextDat3);
                else
                    status = 0;
                end
                % Inkrementiere Zähler
                ndatared = ndatared + 1;
                datacounter = datacounter + 1;
            end % Ende 1. innere Schleife
            idx1(isnan(idx1)) = [];                                        % Lösche dummy Werte
            idx2(isnan(idx2)) = [];
            idx3(isnan(idx3)) = [];
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
                Buf1(1:6,idx1(1:ncols)) = sig;                             % Spannungen
                Buf1(7:12,idx1(1:ncols)) = eps;                            % Dehnungen
                Buf1(13,idx1(1:ncols)) = A(1,:);                           % Durchlaufzähler
                Buf2(1:6,idx2(1:ncols)) = sig;                             % Spannungen
                Buf2(7:12,idx2(1:ncols)) = eps;                            % Dehnungen
                Buf2(13,idx2(1:ncols)) = A(1,:);                           % Durchlaufzähler
                Buf3(1:6,idx3(1:ncols)) = sig;                             % Spannungen
                Buf3(7:12,idx3(1:ncols)) = eps;                            % Dehnungen
                Buf3(13,idx3(1:ncols)) = A(1,:);                           % Durchlaufzähler
            end
            % Buffer ist voll jetzt Zählen
            % 2. innere Schleife über die neuen Werte
            K1 = Nf1(RES1(IZ1));                                           % Zeiger K auf Startindex setzten (Mode I)
            K2 = Nf2(RES2(IZ2));                                           % (Mode II)
            K3 = Nf3(RES3(IZ3));                                           % (Mode III)
            while K1 ~= nextDat1 && ...
                  K2 ~= nextDat2 && ...
                  K3 ~= nextDat3
                % Rainflow
                % ... Mode I
                [IZ1,IR1,RES1,counter1,Nf1,Vg1,P,p] = obj.hcm(...
                    K1,Buf1,Nf1,Vg1,IZ1,IR1,RES1,counter1,nextDat1,7,1,P,p);
                % ... Mode II
                [IZ2,IR2,RES2,counter2,Nf2,Vg2,P,p] = obj.hcm(...
                    K2,Buf2,Nf2,Vg2,IZ2,IR2,RES2,counter2,nextDat2,12,2,P,p);
                % ... Mode III
                [IZ3,IR3,RES3,counter3,Nf3,Vg3,P,p] = obj.hcm(...
                    K3,Buf3,Nf3,Vg3,IZ3,IR3,RES3,counter3,nextDat3,10,3,P,p);
                % Zeiger auf nächsten Wert setzten
                K1 = Nf1(K1);
                K2 = Nf2(K2);
                K3 = Nf3(K3);
            end % Ende 2. innere Schleife
            
            % Freigewordene Stellen (geschlossene Hyst ermitteln)
            empty1 = 1;
            i = nextDat1;
            while Nf1(i) ~= RES1(IR1)
                empty1 = empty1 + 1;
                i = Nf1(i);
            end
            empty2 = 1;
            i = nextDat2;
            while Nf2(i) ~= RES2(IR2)
                empty2 = empty2 + 1;
                i = Nf2(i);
            end
            empty3 = 1;
            i = nextDat3;
            while Nf3(i) ~= RES3(IR3)
                empty3 = empty3 + 1;
                i = Nf3(i);
            end
            
            % Falls der Buffer voll ist aber Daten nicht zu ende sind
            % überlauf anzeigen
            if (empty1 <= 1 || empty2 <= 1 || empty3 <= 1) && 2*nBuf < nBufMax
%               msg = 'Bufferüberlauf in Rainflow PRAM:Buffergröße verdoppeln';
%               warning(msg)
                % verdopple Buffergröße
                Buf1 = [Buf1,zeros(nSp,nBuf)];                             % Neuer Buffer
                Buf2 = [Buf2,zeros(nSp,nBuf)];                             % Neuer Buffer
                Buf3 = [Buf3,zeros(nSp,nBuf)];                             % Neuer Buffer
                Vg1 = [Vg1,nBuf:2*nBuf-1];                                 % Neuer Vorgänger
                Vg1(nBuf + 1) = nextDat1;
                Vg1(Nf1(nextDat1)) = 2*nBuf;   
                Vg2 = [Vg2,nBuf:2*nBuf-1];                                 % Neuer Vorgänger
                Vg2(nBuf + 1) = nextDat2;
                Vg2(Nf2(nextDat2)) = 2*nBuf;
                Vg3 = [Vg3,nBuf:2*nBuf-1];                                 % Neuer Vorgänger
                Vg3(nBuf + 1) = nextDat3;
                Vg3(Nf3(nextDat3)) = 2*nBuf;
                Nf1 = [Nf1,nBuf+2:2*nBuf+1];                               % Neuer Nachfolger 
                Nf1(2*nBuf) = Nf1(nextDat1); 
                Nf1(nextDat1) = nBuf + 1; 
                Nf2 = [Nf2,nBuf+2:2*nBuf+1];                               % Neuer Nachfolger 
                Nf2(2*nBuf) = Nf2(nextDat2); 
                Nf2(nextDat2) = nBuf + 1;
                Nf3 = [Nf3,nBuf+2:2*nBuf+1];                               % Neuer Nachfolger 
                Nf3(2*nBuf) = Nf3(nextDat3); 
                Nf3(nextDat3) = nBuf + 1; 
                empty1 = empty1 + nBuf;
                empty2 = empty2 + nBuf;
                empty3 = empty3 + nBuf;
                nBuf = 2*nBuf;
            elseif empty1 <= 1 || empty2 <= 1 || empty3 <= 1
                msg = 'Bufferüberlauf Programm Rainflow PZ';
                error(msg);
            end
        end % Ende Schleife über alle Werte

        
        %==================================================================
        % Aussortieren leere PZ Werte
        P = P(:,P(1,:) ~= 0);  
        
        %==================================================================
        % clear Memory für Dehnungsvorgeschichte
        obj.exmax = 0;
        obj.exmin = 0;
        obj.exop_alt = 0;
        obj.exop_ein = 0;
        
        %==================================================================
        % Schließe Lastdatei
        fclose(fidLoad);
        
    end % Ende Rainflow
    
    % ... hcm algorithmus
    function [IZ,IR,RES,counter,Nf,Vg,P,pz] = hcm(obj,...
                       K,Data,Nf,Vg,IZ,IR,RES,counter,nextDat,cc,mode,P,pz)
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
                        % ... Display ausgabe Rissöffnungsdehnungen
                        %                     fprintf('%.6f;    %.6f;   ', exop_ein, exop_alt)
                        % ... Schneide Schwingspiel aus
                        [SUBDATA,In,Jn,Kn] = CutOutHyst2(Data,I,J,K,Nf);                        
                        % ... Schädigungsrechnung
                        if mode == 1                           
                            Pz = obj.kurzriss_mode1(SUBDATA,In,Jn,Kn);                           
                        elseif mode == 2                           
                            Pz = obj.kurzriss_mode23( mode,SUBDATA,In,Jn,Kn);                           
                        elseif mode == 3                            
                            Pz = obj.kurzriss_mode23( mode,SUBDATA,In,Jn,Kn);                           
                        end  
                        % ... Speichern mode und Schwingspiel
                        if Pz > obj.PZD0 * 0.1
                            pz = pz + 1;
                            P(:,pz) = [mode;counter;Pz;Data(13,K)];
                        end
                        % ... Rissöffnungsdehnung abklingen lassen !!!
                        obj.exop_alt =  obj.exop_ein + (obj.exop_alt - obj.exop_ein)...
                                                     * exp(-obj.FA/obj.Q*Pz^obj.mJ);
                        % ... Display ausgabe Rissöffnungsdehnungen
%                         fprintf('%.12f;    %.12f\n', obj.exop_ein, obj.exop_alt)
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
    
    % ... Pz Parameter bei Mode I
    function P = kurzriss_mode1(obj,DATA,I,J,K)
        % Funktion zum berechnen des Schädigungsparameters PZ für Mode I
        % Rissöffnung
        %
        %
        % INPUT:
        %
        %  Aus Rainflow/Kerbnäherung:
        %  DATA      -> Spannungen und Dehnungen
        %   I        -> Index Start Hysterese
        %   J        -> Index UKP Hysterese
        %   K        -> Index Ende Hysterese
        %
        %  (Rest sind objekteigenschaften)
        %  Materialparameter:
        %  sigG      -> Fließspannung
        %
        %  Geometrie:
        %  Y1A0      -> Geometrie Faktor Riss
        %
        %  Mittelspannungsempfindlichekeit
        %  M_sig     -> Mittelspannungsempfindlichekeit
        %
        %  Geschichtsvariablen:
        %  exop_alt  -> Rissöffnungsdehnung aus Vorgeschichte, wird hier
        %               aktuallisiert
        %  exmax     -> Bisher aufgetretene Maximaldehnung, wird hier
        %               aktuallisiert
        %  exmin     -> Bisher aufgetretene Minimaldehnung, wird hier
        %               aktuallisiert
        %
        % OUTPUT:
        %
        %    P      -> Schädigungsparameter
        %   ...     -> aktuallisierte Eigenschaften der 
        %              Rissöffnungsdehnungen 
        %__________________________________________________________________
        
        
        % -----------------------------------------------------------------
        % Maximale und minimale Normalspannung
        [sx_max,imax] = max(DATA(1,I:K));
        imax = imax + I - 1;
        
        sx_min = min(DATA(1,I:K));
               
        % -----------------------------------------------------------------
        % obere von Mises Vergleichsspannung und obere wirksame 
        % Vergleichsspannung beide mit dem Vorzeichen der größten
        % Normalspannung
        [sv_vz, sw_vz] = obj.mises_vorzeichen(DATA,imax);
        
        
        % -----------------------------------------------------------------
        % Rissöffnungsspannung bestimmen
        [sxop] = obj.newman_max(sx_max,sx_min,sw_vz,sv_vz,obj.sigF,obj.M_sig,obj.ators);
        
        
        % -------------------------------------------------------------------------
        % Indizes der Hysterese festlegen (unterscheide Auf- und Absteigende Äste)
        % und stehender/hängender Hysteresen
        if DATA(1,I) < DATA(1,J) % Start der Hysterese kleine als UKP
            % ... stehende Hysterese
            iu_auf = I;
            io_auf = J;
            io_ab = J;
            iu_ab = K;
        else
            % ... hängende Hysterese
            io_ab = I;
            iu_ab = J;
            iu_auf = J;
            io_auf = K;
        end
        
        
        % -------------------------------------------------------------------------
        % Rissöffnungspunkt und -dehnung am aufsteigenden Ast suchen
        obj.exop_ein = obj.epsopen(DATA,iu_auf,io_auf,sxop,sx_max);
        
        
        % -------------------------------------------------------------------------
        % Reihenfolgeeffekte
        
        exo = DATA(7,io_ab);            % Dehnung oben
        exu = DATA(7,iu_auf);           % Dehnung unten
        sxa = 0.5 * (sx_max - sx_min);  % Spannungsamplitude
        % ... abfangen gedrehte Hysteresen
%         if exu > exo
%             dummy = exu;
%             exu = exo;
%             exo = dummy;
%         end
        
        if exo <= obj.exop_alt
            % ... Riss komplett zu
            exop = exo;
            
        elseif exo > obj.exmax || exu < obj.exmin % Über/Unterschreitung der Maximalwerte
            % !!!! HIER SIND MÖGLICHERWEISE TOLERANZEN HILFREICH UM NUMERISCHE
            % UNTERSCHIEDE ABZUFANGEN !!!!!!
            % ... neu Setzten der Extrema
            obj.exmax = exo;
            obj.exmin = exu;
            exop = obj.exop_ein;
            obj.exop_alt = obj.exop_ein;
            
        elseif obj.exop_ein >= obj.exop_alt
            % ... Schwingspiel anhand der aktuellen(alten) Öffnungsdehnung bewerten und
            % anschließend neue nehmen
            exop = obj.exop_alt;
            
        elseif sxa >= 0.4*obj.sigF
            % ... großes Schwingspiel setzt Öffnungsdehnung auf niedrigeres Niveau
            exop = obj.exop_ein;
            obj.exop_alt = obj.exop_ein;
            
        else
            % ... Rissöffnungsdehnung aus Vorgeschichte
            exop = obj.exop_alt;
            
        end
        
        % ... Display ausgabe verwendete Öffnungsdehnung
        % fprintf('%.6f;   ', exop)
        
        % -------------------------------------------------------------------------
        % Rissschließdehnung
        
        excl = exop;
        
        % -------------------------------------------------------------------------
        % Rissschließpunkt und -spannung am absteigenden Ast suchen
        
        [sxcl,icl] = obj.sigclose(DATA,iu_ab,io_ab,excl,sx_max,sx_min);
        
        % ... plote Hysterese mit Schließpunkt
        % plot(DATA(7,:),DATA(1,:),'b','LineWidth',1.4)
        % hold on
        % plot(excl,sxcl,'.k','MarkerSize',10);
        % pause(5)
        % -------------------------------------------------------------------------
        % Integration von DWxeff am absteigenden Ast
        
        Wxeff = obj.verzerrungsenergie(DATA,1,7,io_ab,icl);
        
        % Korrektur der Verzerrungsenergie, da über den exakten (aus diskreten
        % punkten interpolierten) Schließpunkt drüber integriert wurde
        Wxeff = Wxeff - 0.5 * ( DATA(1,icl)+sxcl-2*DATA(1,io_ab) ) * ...
            ( DATA(7,icl) - excl );
        
        % Ausschlißen Wx < 0
        Wxeff = max([0,Wxeff]);
        
        % -------------------------------------------------------------------------
        % PZ-Wert bestimmen
        
        P = 2*pi*Wxeff*obj.Y1A0^2;
        % P = 6.2632*Wxeff*Y1A0^2;
        
    end % Ende Kurzrissmode I
    
    % ... Pz Parameter bei Mode II & III
    function P = kurzriss_mode23(obj,mode,DATA,I,J,K)
        % Funktion zum berechnen des Schädigungsparameters PZ für Mode II 
        % und III Rissöffnung
        %
        % INPUT:
        %
        %  Aus Rainflow/Kerbnäherung:
        %  mode      -> Rissöffnungsmode (2) oder (3)
        %  DATA      -> Spannungen und Dehnungen
        %   I        -> Index Start Hysterese
        %   J        -> Index UKP Hysterese
        %   K        -> Index Ende Hysterese
        %  (Rest sind Eigenschaften der Klasse)
        %  Materialparameter:
        %  sf        -> Fließspannung
        %  takt      -> Schubspannung für Reibkontakt der Rissufer
        %  mu        -> Faktor zur Berücksichtigung des Einfluss der Normalspannung
        %               auf Rissuferkontakt
        %  nu        -> Querdehnzahl
        %
        %  Geometrie:
        %  YA0       -> Geometrie Faktor Riss
        %              (Y2A0 für Mode II und Y3A0 für Mode III)
        %
        % OUTPUT:
        %
        %    PZ      -> Schädigungsparameter
        %
        %__________________________________________________________________________
        
        % -------------------------------------------------------------------------
        % Festlegen der auszuwertenden Größen je nach Mode am Punkt A !!!
        if mode == 2
            zs = 6;         % Zeile der sxz Spannungen
            ze = 12;        % Zeile der gxz Dehnungen
        elseif mode == 3
            zs = 4;         % Zeile der sxy Spannungen
            ze = 10;        % Zeile der gxy Dehnungen
        else
            msg = 'Fehler in Kurzrissmodell Mode 23: kein gültiger Mode angegeben';
            error(msg)
        end
        
        
        % -------------------------------------------------------------------------
        % Maximale und minimale Schubspannung
        % maximal
        t_max = max(DATA(zs,I:K));
        % minimal
        t_min = min(DATA(zs,I:K));
        % Schwingweite Schubspannungen
        dtau = abs(t_max - t_min);
        % Schwingweite der Schubverzerrungen
        dgam1 = abs(DATA(ze,I) - DATA(ze,J));
        dgam2 = abs(DATA(ze,J) - DATA(ze,K));
        
        % -------------------------------------------------------------------------
        % Mittelwert der Normalspannungen (gewichtet mit Dehnungsinkrementen)
        
        snm = 0.5 * sum( ( DATA(1,I+1:K) + DATA(1,I:K-1) ) .* ...
            abs( DATA(ze,I+1:K) - DATA(ze,I:K-1) ) );
        snm = snm/(dgam1+dgam2);
        
        
        % -------------------------------------------------------------------------
        % Ueff aus Mittenrissscheibe mit Dugdale Modell
        
        % Schubfließspannung 
%         tf = obj.sigF/sqrt(3);
        tf = obj.tauF;
        
        % Reibspannung
        tr = max([0,obj.takt-obj.mu*snm]);
        
        % Unterscheide geöffnete/geschlossene Riss
        if dtau <= 2 * tr || dtau < 1e-6
            % ... Riss komplett zu
            Ueff = 0;
            
        elseif dtau > 2 * tf || tr <= 0
            % ... Riss komplett auf
            Ueff = 1;
        else
            % ... Riss teilweise auf (Ueff mit Mode 2 Dugdale)
            Ueff = log( cos( pi/2 * (dtau/2-tr)/(tf-tr) ) )/...
                log( cos( pi/2 * dtau/(2*tf)) );
            
            if Ueff > 1
                Ueff = 1;
            elseif Ueff < 0
                Ueff = 0;
            end
            
        end
        
        % -------------------------------------------------------------------------
        % Effektive Verzerrungsenergie
        
        % auf aufsteigendem und absteigendem Ast und dann Mittelwert bilden
        W1 = obj.verzerrungsenergie(DATA,zs,ze,I,J);
        W1 = max([0,W1]);
        W2 = obj.verzerrungsenergie(DATA,zs,ze,J,K);
        W2 = max([0,W2]);
        
        % Effektive Verzerrungsenergie
        Weff = 0.5 * ( W1 + W2 ) * Ueff;
        
        
        % -------------------------------------------------------------------------
        % PZ bestimmen
        
        % unterscheide Moden
        if mode == 2
            % ... Mode II am Punkt A
            P = pi/(1+obj.nu) * obj.Y2A0^2 * Weff;
        else
            % ... Mode III am Punkt A
            P = pi * Weff * obj.Y3A0^2;
        end
    end % Ende Kurzrissmode II & III
    
    % ... Lebensdauerrechnung
    function [DL,SSP] = lebensdauer(obj,P,varargin)
        % Funktion rechnet Lebensdauern aus Schädigungsparametern
        %
        % INPUT:
        % P                 - Schädingungsparameter
        %                      (1.Zeile) -> Rissöffnungsmode
        %                      (2.Zeile) -> gezählte Schwingspiele
        %                      (3.Zeile) -> PZ Parameter
        %                      (4.Zeile) -> Durchlauf in dem Ssp gezählt wurde
        % varargin          - Variabler input für output datei
        %
        % OUTPUT:
        % DL             - Durchläufe
        % SSP            - Schwingspiele
        %__________________________________________________________________
        % -----------------------------------------------------------------
        % Durchläufe
        ndl = ceil(max(P(4,:)));              % Anzahl Durchläufe
        
        % -----------------------------------------------------------------
        % Anfangswerte Riss
        ai = obj.a0;
        ci = obj.c0;
        da = 0.1*ai;     % Inkrementelle Risslänge mit konstantem a/c
        Dsum = 0;        % Schadenssumme
        DsumModes = zeros(3,1); % Schadenssumme der einzelnen Moden
        azc = ai/ci;     % Verhältniss Risslängen a/c
        
        % -------------------------------------------------------------------------
        % Organisiere Output
        if nargin == 3
            % ... output erstellen
            outopt = 1;
            % ... Speicherzähler 
            speicherzahler = 1;
            % ... Speicher vorbehalten
            obj.speicherrifo = zeros(1000,5);
            % ... Speicher Startwerte
            obj.speicherrifo(speicherzahler,:) = [0,0,ai*1000,ci*1000,ai/ci];
        else
            % ... kein output erstellen
            outopt = 0;
        end
        
        % -------------------------------------------------------------------------
        % PZ-Wöherlinie
        [PZD,ND_pz,dc] = obj.ink_rifo(ai,ci,da);
        
        % -------------------------------------------------------------------------
        % Abfangen leeres PZ
        if isempty(P)
            DL = obj.dauerfest;                       % Dummywert für dauerfest
            SSP = obj.dauerfest;
%             msg = 'PZ:keine Schädigung identifieziert';
%             warning(msg)
            return;
        end
        
        % -------------------------------------------------------------------------
        % Alle Durchläufe durchgehen (außer den Letzten)
        i = 1;
        while ci < obj.ce && P(4,i) < ndl - 1
            
            % ... Variablen
            mode = P(1,i);                     % Rissöffnungsmode
            pz = P(3,i);                       % aktueller PZ-Wert
            
            % ... lineare Schadensakkumulation
            [Dsum,Dakt] = obj.DamAkk(pz,obj.mJ,ND_pz(mode),PZD(mode),Dsum);
            
            % ... update RissLänge
            ai = ai + da*Dakt;
            ci = ci + dc(mode)*Dakt;
            DsumModes(mode) = DsumModes(mode) + Dakt;
            
            % ... update PZ-WL wenn Dsum = 1
            if Dsum >= 1
                % ... Schadenssumme wieder auf Null
                Dsum = 0;
                DsumModes = zeros(3,1);
                              
                % ... Update PZ-WL
                [PZD,ND_pz,dc,azc] = obj.ink_rifo(ai,ci,da);
                
                % ... Speichern Werte
                if outopt
                    % ... Speicher Startwerte
                    speicherzahler = speicherzahler + 1;
                    obj.speicherrifo(speicherzahler,:) = [P(4,i),P(2,i),ai*1000,ci*1000,ai/ci];               
                end
            end
            
            % ... inkrement i
            i = i + 1;
            
        end % Ende Schleife über alle PZ Werte
        
        ilast = i;                     % Zeiger auf start des letzten durchlaufs
        c1 = ci;                       % Merken Risslänge
        
        % -------------------------------------------------------------------------
        % Schauen obs kaputt is
        if ci >= obj.ce
            % ... kaputt
            DL = P(4,i-1);                  % Durchläufe bis Kaputt
            SSP = i - 1;                    % Gesamt Schwingspiele (whs sinnlose info)
            
        elseif ci <= obj.c0 && ndl > 1
            % ... dauerfest (vllt kommt noch was im letzten Durchlauf, aber whs net)
            DL = obj.dauerfest;                       % Dummywert für dauerfest
            SSP = obj.dauerfest;
        else
            % ... letzten Durchlauf aufbringen
            while ci < obj.ce && i <= size(P,2)
                
                % ... Variablen
                mode = P(1,i);                     % Rissöffnungsmode
                pz = P(3,i);                       % aktueller PZ-Wert
                
                % ... lineare Schadensakkumulation
                [Dsum,Dakt] = obj.DamAkk(pz,obj.mJ,ND_pz(mode),PZD(mode),Dsum);
                
                % ... update RissLänge
                ai = ai + da*Dakt;
                ci = ci + dc(mode)*Dakt;
                DsumModes(mode) = DsumModes(mode) + Dakt;
                
                % ... update PZ-WL wenn Dsum = 1
                if Dsum >= 1
                    % ... Schadenssumme wieder auf Null
                    Dsum = 0;
                    DsumModes = zeros(3,1);
                    
                    % ... Update PZ-WL
                    [PZD,ND_pz,dc,azc] = obj.ink_rifo(ai,ci,da);
                    
                    % ... Speichern Werte
                    if outopt
                        % ... Speicher Startwerte
                        speicherzahler = speicherzahler + 1;
                        obj.speicherrifo(speicherzahler,:) = [P(4,i),P(2,i),ai*1000,ci*1000,ai/ci];
                    end
                end
                
                % ... inkrement i
                i = i + 1;
                
            end % Ende Schleife über alle PZ Werte
            
            % ... Schauen obs kaputt is
            if ci >= obj.ce
                % ... Kaputt
                DL = P(4,i-1);                  % Durchläufe bis Kaputt
                SSP = i - 1;                    % Gesamt Schwingspiele (whs sinnlose info)
            elseif ci <= c1
                % ... dauerfest
                DL = obj.dauerfest;                       % Dummywert für dauerfest
                SSP = obj.dauerfest;
            else
                % letzten Durchlauf solange aufbringen bis kaputt
                % DL init & SSP init
                DL = ndl;
                SSP = i;
                % Zeiger zurücksetzten
                i = ilast;
                % merke Schadenssumme  
                rifoDL = 1;            % schaut ob in einem DL Dsum = 1 erreicht wird
                Dsum1 = Dsum;
                DsumModes1 = DsumModes;
%                 if ilast > 1
%                     dDL = P(4,end)-P(4,ilast-1);
%                 else
%                     dDL = P(4,end)-P(4,ilast);
%                 end
                dDL = 1;
                % Schleife
                while ci < obj.ce && DL < obj.dauerfest
                    
                    % ... Variablen
                    mode = P(1,i);                     % Rissöffnungsmode
                    pz = P(3,i);                       % aktueller PZ-Wert
                    
                    % ... lineare Schadensakkumulation
                    [Dsum,Dakt] = obj.DamAkk(pz,obj.mJ,ND_pz(mode),PZD(mode),Dsum);
                    
                    % ... update RissLänge
                    ai = ai + da*Dakt;
                    ci = ci + dc(mode)*Dakt;
                    DsumModes(mode) = DsumModes(mode) + Dakt;
                    
                    % update Schwingspiele
                    SSP = SSP + 1;
                    if i > 1
                        DL = DL + P(4,i) - P(4,i-1);
                    else
                        DL = DL + P(4,i);
                    end
                    
                    % ... update PZ-WL wenn Dsum = 1
                    if Dsum >= 1
                        % ... Schadenssumme wieder auf Null
                        Dsum = 0;
                        DsumModes = zeros(3,1);
                        % ... Update PZ-WL
                        [PZD,ND_pz,dc,azc] = obj.ink_rifo(ai,ci,da);
                        rifoDL = 0;
                        
                        % ... Schreibe Output
                        if outopt
                            % ... Speicher Startwerte
                            speicherzahler = speicherzahler + 1;
                            obj.speicherrifo(speicherzahler,:) = [DL,SSP,ai*1000,ci*1000,ai/ci];
                        end
                    end
                    
                    % ... inkrement i
                    i = i + 1;
                    
                    % ... zurücksetzten i wenn alle Werte des Durchlaufs bearbeitet
                    if i > size(P,2)
                        i = ilast;
                        % ... Überspringe Durchläufe wenn riss nicht soweit
                        % wächst um WL zu updaten
                        if rifoDL && ci < obj.ce
                            dDsum = Dsum - Dsum1;
                            dDsumModes = DsumModes - DsumModes1;
                            ... Schauen ob Dauerfest
                            if dDsum == 0
                                DL = obj.dauerfest;
                            % ... Überspringe Durchläufe
                            else
                                % ... Durchläufe bevor WL geupdatentet wird
                                % oder ci + dc > ce (check hier ai + da > ae)
                                aend = azc * obj.ce;
                                if ai + da > aend
                                    x = floor((aend - ai)/(da * dDsum));
                                else
                                    x = floor((1-Dsum)/dDsum);
                                end    
                                % ... abfangen numerische Fehler
                                x = max([x 0]);
                                % ... update Durchläufe
                                DL = DL + x*dDL;
                                SSP = SSP + x;
                                % ... update Risslänge
                                ai = ai + da * dDsum * x;
                                ci = ci + dc * dDsumModes * x;
                                Dsum = Dsum + dDsum * x;
                                DsumModes = DsumModes + dDsumModes * x;
                            end
                        end
                        % ... Merken aktuelle Werte                        
                        Dsum1 = Dsum;  
                        DsumModes1 = DsumModes;
                        rifoDL = 1;
                    end                   
                end
            end
            
        end
                
    end % Ende Lebensdauerrechnung
    
    % ... Inkrementeller Rissfortschritt
    function [PZD,ND_pz,dc,azc] = ink_rifo(obj,ai,ci,da)
        % Funktion berechnet neue PZ-WL bei ink. Rissfortschrtt
        % Rechnung Vektorisiert unter der Annahme a/c = konst
        %
        % INPUT:
        % ai,ci           - aktuelle Risslängen
        % da              - Inkrementeller Risslängnzuwachs Punkt A
        % (Rest aus Eigenschaften der Klasse)
        % Gsig,Gtau       - bezogene Spannungsgradienten senkrecht zur Oberfläche
        % E,nu            - Elastizitätskonstanten
        % CJ              - Konstante Rissfortschrittsgesetzt (Paris)
        % ZeffthLR        - Langrissschwellenwert für Zeff
        % mJ              - Steigung der PZ-WL (und Exponent Paris Gl.)
        % a0,c0           - mikrostr. Anfangsrisslängen
        % PZD0            - Kurzrissinitiierungsschwellwert (PZ Dauerfestgkeit
        %                   Werkstoffprobe)
        % lstern          - mikrostrukturelle Hilfsgröße
        % Y1A0,Y2A0,Y3A0  - Geometriefunktionen kurzer Riss in homogenem
        %                   Spannungsfeld und a/c =0.9
        %
        % OUTPUT:
        % PZD             - Dauerfestigkeit der einzelnen Moden e R^(1x3)
        % ND              - Nötige Ssp für Rifo um da bei PZD der einzelnen Moden e R^(1x3)
        % dc              - zugehöriger Risslängenzuwachs der Moden am Punkt B e R^(1x3)
        %__________________________________________________________________
        
        
        
        % -------------------------------------------------------------------------
        % a/c Verhältniss
        % ... unsinnige Werte Abfangen
        if ai <= 0
            ai = obj.a0;
        end
        if ci <= 0
            ci = obj.c0;
        end
        % ... a/c
        azc = ai/ci;
        
        % -------------------------------------------------------------------------
        % Risslängenunabhängige Dauerfestigkeit (erstmal ungekerbt)
        PZD = obj.PZD0 * (obj.a0+obj.lstern)/(ai+obj.lstern);
        PZD = [PZD,PZD,PZD]; % erstmal gleich für alle moden
        
        % -------------------------------------------------------------------------
        % a/rho
        azr1 = 0.5 * obj.Gsig*ai;      % aus Gsig = 2/rho_s
        azr2 = obj.Gtau * ai;          % aus Gtau = 1/rho_t
        
        
        % -------------------------------------------------------------------------
        % Geometriefunktionen
%         [Y1A,Y1B] = obj.Y_OFR_TYPA_ZD(azc,azr1);
        Y1A = obj.YIA(azc,azr1);
        Y1B = obj.YIB(azc,azr1);
%         [~,Y2B,Y3A,~] = obj.Y_OFR_TYPA_SH(azc,azr2);
        Y3A = obj.YIIIA(azc,azr2);
        Y2B = obj.YIIB(azc,azr2);
%         [~,~,Y2A,~,~,Y3B] = obj.Y_OFR_TYPB_ZD(azc,azr1);
        Y2A = obj.YIIA(azc,azr1);
        Y3B = obj.YIIIB(azc,azr1);
        
        % -------------------------------------------------------------------------
        % Risslängenunabhängige Dauerfestigkeit (gekerbt)
        % ... Mode I
        PZD(1) = PZD(1) * (obj.Y1A0/Y1A)^2;
        % ... Mode II
        PZD(2) = PZD(2) * (obj.Y2A0/(2*Y2A))^2; % mal 2 wegen hinweis Hertl
        % ... Mode III
        PZD(3) = PZD(3) * (obj.Y3A0/Y3A)^2;
               
        % -------------------------------------------------------------------------
        % Risslängenzuwachs Punkt B (Bei Annahme a/c = konst)
        dc(1) = da * ( Y1B^2/Y1A^2 )^(obj.mJ);
        dc(2) = da * ( Y3B^2/Y2A^2 )^(obj.mJ)*(1+obj.nu)^obj.mJ;     % ! Mode II Punkt A = Mode III Punkt B
        dc(3) = da * ( Y2B^2/Y3A^2 )^(obj.mJ)*(1/(1+obj.nu))^obj.mJ; % ! Mode III Punkt A = Mode II Punkt B
        
        % -------------------------------------------------------------------------
        % Nötige Schwingspiele bis Risswachstum um da bei PZD
        % ... Parameter Integrationsschrittweite
        n = 101;          % Anzahl an Integrationspunkten -> n-1 Inkremente
        h = da/(n-1);     % Schrittweite der Integration
        
        % ... Risslängen
        a = ai:h:ai+da;
        
        % ... a/rho Verhältniss
        azr1 = 0.5 * obj.Gsig * a;
        azr2 = obj.Gtau * a;
        azc = repelem(azc,n);
        
        % ... Geometriefunktionen
%         [Y1A,~] = obj.Y_OFR_TYPA_ZD(azc,azr1);
        Y1A = obj.YIA(azc,azr1);
%         [~,~,Y3A,~] = obj.Y_OFR_TYPA_SH(azc,azr2);
        Y3A = obj.YIIIA(azc,azr2);
%         [~,~,Y2A,~,~,~] = obj.Y_OFR_TYPB_ZD(azc,azr1);
        Y2A = obj.YIIA(azc,azr1);
        
        % ... Funktion F (nicht Def. aus Diss sondern Integrant für ND)
        f1 = (a .* (Y1A./obj.Y1A0).^2 ).^(-obj.mJ);
        f2 = (a .* (2*Y2A./obj.Y2A0).^2 ).^(-obj.mJ); % mal 2 wegen Hinweis Hertl
        f3 = (a .* (Y3A./obj.Y3A0).^2 ).^(-obj.mJ);
        
        % -------------------------------------------------------------------------
        % Integration von F
        F1 = obj.simpson(f1,h,n);
        F2 = obj.simpson(f2,h,n);
        F3 = obj.simpson(f3,h,n);
        
        % -------------------------------------------------------------------------
        % Anzahl Ssp bei PZD
        ND_pz = zeros(1,3);
        ND_pz(1) = F1/(obj.CJ*PZD(1)^obj.mJ);
        ND_pz(2) = F2/(obj.CJ*PZD(2)^obj.mJ);
        ND_pz(3) = F3/(obj.CJ*PZD(3)^obj.mJ);
    end % Ende Inkrementeller Rissfortschritt
    
    % ... Plote PWL
    function plotPWL(obj,varargin)
        % Plotet die PZ-Wöhlerlinie
        % PZ = (N/Q)^-(1/mJ)
        % varargin(1) - axis handle
        % varargin    - beliebiege Plot Optionen
        %------------------------------------------------------------------
        % varargin = ax
        if nargin == 1
            figure, grid on, hold on
            ax = gca;
        else
            ax = varargin{1};
        end
        % ... Lebensdauerpunkte [1,Dauerfestigkeit,irgendwo über Dauerfest]
        N = [1, obj.ND0, obj.ND0*1000];
        % ... PZ Werte
        P = [(N(1)/obj.Q)^-(1/obj.mJ),obj.PZD0,obj.PZD0];
        % ... Plot
        plot(ax,N,P,varargin{2:end});
        set(gca,'XScale','log','YScale','log');
        xlabel('N'), ylabel('P')
        grid on
    end
    
    % ... Geometriefaktoren
    function [YIIA,YIIB,YIIIA,YIIIB] = Y_OFR_TYPA_SH(obj)
        % Funktion gibt risslängenabhängige Geometriefunktionen in abhängigkeit des
        % Risslängenverhältnisses a/c und Geometrieverhältniss a/rho
        %
        % 
        % Geometriefunktionen sind für Oberflächenriss unter Schubbelastung
        % Risstiefe a=1, Schubspannung txy=1
        % Punkt A - tiefster Punkt
        % Punkt B - bei theta = 10°
        %
        % Werte werden interpoliert aus Stützstellen berechnet von Hertel.
        % Siehe Diss Hertl.
        %
        % INPUT:
        % azc        - Verhältniss a zu c (Stützstellen 1,3/4,2/3,1/2,1/4)
        % azr        - Verhältniss a zu rho (fiktiver Kerbradius) (Stützstellen 0,1/4,1/2,1,2,4)
        %
        % OUTPUT:
        %  Y2A,Y2B   - Geometriefunktionen Mode II Stelle A und B 
        %  Y3A,Y3B   - Geometriefunktionen Mode III Stelle A und B
        %
        %__________________________________________________________________________
        
        % -------------------------------------------------------------------------
        % Stützstellen (aus Diss Hertel)
        % ... a zu c
        ac = [ 0.25 0.5 0.667 0.75 1 ];
        % ... a zu rho
        ar = [ 0 0.25 0.5 1 2 4];
        
        % ... Geometriefunktionen aus Abaqus
        %  a/c= 0.25    0.5     0.667    0.75    1            a/rho
        f3a = [0.8982, 0.7588, 0.6753, 0.6370, 0.5376;...   %  0
            0.7114, 0.5997, 0.5332, 0.5030, 0.4238;...   %  0.25
            0.6105, 0.5150, 0.4578, 0.4319, 0.3641;...   %  0.5
            0.5014, 0.4232, 0.3764, 0.3551, 0.2994;...   %  1
            0.3965, 0.3350, 0.2980, 0.2812, 0.2373;...   %  2
            0.3042, 0.2574, 0.2291, 0.2162, 0.1825];...  %  4
            
        f2b = [0.5129, 0.6562, 0.6883, 0.6956, 0.6973;...   %  0
            0.4191, 0.5339, 0.5604, 0.5659, 0.5681;...   %  2.5
            0.3628, 0.4611, 0.4840, 0.4890, 0.4912;...   %  0.5
            0.3008, 0.3808, 0.3997, 0.4039, 0.4060;...   %  1
            0.2403, 0.3027, 0.3178, 0.3211, 0.3230;...   %  2
            0.1863, 0.2334, 0.2449, 0.2475, 0.2491];     %  4
        
        % f3b = [0.3405, 0.2808, 0.2533, 0.2411, 0.2077;...   %  0
        %        0.2777, 0.2277, 0.2051, 0.1951, 0.1681;...   %  0.25
        %        0.2403, 0.1965, 0.1769, 0.1684, 0.1451;...   %  0.5
        %        0.1991, 0.1621, 0.1460, 0.1389, 0.1198;...   %  1
        %        0.1589, 0.1287, 0.1159, 0.1104, 0.0952;...   %  2
        %        0.1232, 0.0992, 0.0893, 0.0850, 0.0734];     %  4
        
        % -------------------------------------------------------------------------
        % erstelle Netz
%         [AC,AR] = meshgrid(ac,ar);
        [AC,AR] = ndgrid(ac,ar);
        % -------------------------------------------------------------------------
        % Interpolation von zwischenwerten (nur Werte die auch gebraucht werden,
        % wegen zeit)
        
        % ... Modus II
        YIIA = 0;
        YIIB = griddedInterpolant(AC, AR, f2b',obj.InterpMethod);
        % ... Modus III
        YIIIA = griddedInterpolant(AC, AR, f3a',obj.InterpMethod);
        YIIIB = 0; % Eigentlich nicht Null aber Wert wird nicht gebraucht
%         YIIIB = griddedInterpolant(AR, AC, f3b,obj.InterpMethod);
        
        
    end % Ende Funktion
    function [YIA,YIB] = Y_OFR_TYPA_ZD(obj)
        % Funktion gibt risslängenabhängige Geometriefunktionen in abhängigkeit des
        % Risslängenverhältnisses a/c und Geometrieverhältniss a/rho
        %
        % Geometriefunktionen sind für Oberflächenriss unter Zugbeanspruchung
        % Risstiefe a=1, Zugspannung sx=1
        % Punkt A - tiefster Punkt
        % Punkt B - bei theta = 10°
        %
        % Werte werden interpoliert aus Stützstellen berechnet von Hertel.
        % Siehe Diss Hertl.
        %
        % INPUT:
        % azc        - Verhältniss a zu c (Stützstellen 1,3/4,2/3,1/2,1/4)
        % azr        - Verhältniss a zu rho (fiktiver Kerbradius) (Stützstellen 0,1/4,1/2,1,2,4)
        %
        % OUTPUT:
        %  Y1A,Y1B   - Geometriefunktionen Mode I Stelle A und B
        %
        %__________________________________________________________________________
        
        
        % -------------------------------------------------------------------------
        % Stützstellen (aus Diss Hertel)
        % % ... a zu c
        ac = [ 0.25 0.5 0.667 0.75 1 ];
        % ... a zu rho
        ar = [ 0 0.25 0.5 1 2 4];
        
        % ... Geometriefunktionen aus Abaqus
        %  a/c= 0.25    0.5     0.667    0.75    1            a/rho
        f1a = [1.0270, 0.8846, 0.7995, 0.7605, 0.6601;...   %  0
            0.6629, 0.5704, 0.5143, 0.4889, 0.4256;...   %  0.25
            0.5223, 0.4497, 0.4054, 0.3853, 0.3359;...   %  0.5
            0.3947, 0.3402, 0.3065, 0.2913, 0.2543;...   %  1
            0.2900, 0.2503, 0.2254, 0.2142, 0.1871;...   %  2
            0.2095, 0.1810, 0.1629, 0.1548, 0.1353];     %  4
        f1b = [0.5789, 0.6883, 0.7165, 0.7227, 0.7244;...   %  0
            0.3867, 0.4596, 0.4782, 0.4823, 0.4857;...   %  0.25
            0.3053, 0.3635, 0.3782, 0.3816, 0.3848;...   %  0.5
            0.2311, 0.2753, 0.2866, 0.2891, 0.2920;...   %  1
            0.1704, 0.2024, 0.2109, 0.2128, 0.2152;...   %  2
            0.1236, 0.1462, 0.1524, 0.1538, 0.1557];     %  4
        
        % -------------------------------------------------------------------------
        % erstelle Netz
%         [AC,AR] = meshgrid(ac,ar);
        [AC,AR] = ndgrid(ac,ar);
        
        % -------------------------------------------------------------------------
        % Interpolation von zwischenwerten
        
        YIA = griddedInterpolant(AC, AR, f1a', obj.InterpMethod);
        YIB = griddedInterpolant(AC, AR, f1b', obj.InterpMethod);
    end % Ende Funktion
    function [YIA,YIB,YIIA,YIIB,YIIIA,YIIIB] = Y_OFR_TYPB_ZD(obj)
        % Funktion gibt risslängenabhängige Geometriefunktionen in abhängigkeit des
        % Risslängenverhältnisses a/c und Geometrieverhältniss a/rho
        %
        % Geometriefunktionen sind für Oberflächenriss der um 45° geneigt ist
        % unter Zugbeanspruchung
        % Risstiefe a=1, Zugspannung sx=1
        % Punkt A - tiefster Punkt
        % Punkt B - bei theta = 10°
        %
        % Werte werden interpoliert aus Stützstellen berechnet von Hertel.
        % Siehe Diss Hertl.
        %
        % INPUT:
        % azc        - Verhältniss a zu c (Stützstellen 1,3/4,2/3,1/2,1/4)
        % azr        - Verhältniss a zu rho (fiktiver Kerbradius) (Stützstellen 0,1/4,1/2,1,2,4)
        %
        % OUTPUT:
        %  Y1A,Y1B   - Geometriefunktionen Mode I Stelle A und B
        %  Y2A,Y2B   - Geometriefunktionen Mode II Stelle A und B
        %  Y3A,Y3B   - Geometriefunktionen Mode III Stelle A und B
        %
        %__________________________________________________________________________
        
        % -------------------------------------------------------------------------
        % Stützstellen (aus Diss Hertel)
        % ... a zu c
        ac = [ 0.25 0.5 0.667 0.75 1 ];
        % ... a zu rho
        ar = [ 0 0.25 0.5 1 2 4];
        
        % ... Geometriefunktionen aus Abaqus
        %  a/c= 0.25    0.5     0.667    0.75    1            a/rho
        % f1a = [0.6533,	0.5715,	0.5119,	0.4830,	0.4059;...   %  0
        %        0.4653,	0.4089,	0.3658,	0.3449,	0.2892;...   %  0.25
        %        0.3746,	0.3311,	0.2962,	0.2792,	0.2340;...   %  0.5
        %        0.2865,	0.2555,	0.2286,	0.2154,	0.1804;...   %  1
        %        0.2114,	0.1906,	0.1704,	0.1606,	0.1343;...   %  2
        %        0.1528,	0.1392,	0.1244,	0.1171,	0.0978];...  %  4
        
        f2a = [0.3539,	0.3312,	0.3214,	0.3173,	0.3049;...   %  0
            0.2511,	0.2331,	0.2261,	0.2233,	0.2146;...   %  0.25
            0.2045,	0.1884,	0.1826,	0.1803,	0.1734;...   %  0.5
            0.1593,	0.1453,	0.1405,	0.1387,	0.1334;...   %  1
            0.1203,	0.1084,	0.1045,	0.1031,	0.0992;...   %  2
            0.0889,	0.0791,	0.0760,	0.0749,	0.0721];...  %  4
            
        % f1b = [0.3286, 0.3713, 0.3899, 0.3959, 0.4034;...   %  0
        %        0.2424, 0.2737, 0.2875, 0.2920, 0.2976;...   %  0.25
        %        0.1973, 0.2223, 0.2337, 0.2374, 0.2419;...   %  0.5
        %        0.1531, 0.1718, 0.1806, 0.1835, 0.1871;...   %  1
        %        0.1152, 0.1281, 0.1347, 0.1368, 0.1396;...   %  2
        %        0.0850, 0.0934, 0.0981, 0.0998, 0.1018];     %  4
        %
        % f2b = [0.06426, -0.01201 , -0.03301 , -0.03943 , -0.04914;...   %  0
        %        0.04761, -0.008294, -0.02379 , -0.02854 , -0.03580;...   %  0.25
        %        0.03855, -0.006680, -0.01926 , -0.02314 , -0.02906;...   %  0.5
        %        0.02954, -0.005167, -0.01486 , -0.01785 , -0.02244;...   %  1
        %        0.02179, -0.003902, -0.01108 , -0.01330 , -0.01673;...   %  2
        %        0.01575, -0.002904, -0.008079, -0.009693, -0.01219];     %  4
        
        f3b = [-0.1747, -0.2336, -0.2513, -0.2568, -0.2652;...   %  0
            -0.1287, -0.1717, -0.1847, -0.1887, -0.1946;...   %  0.25
            -0.1044, -0.1395, -0.1501, -0.1533, -0.1581;...   %  0.5
            -0.0805, -0.1077, -0.1159, -0.1185, -0.1222;...   %  1
            -0.0599, -0.0802, -0.0864, -0.0884, -0.0912;...   %  2
            -0.0436, -0.0583, -0.0630, -0.0644, -0.0665];     %  4
        
        % -------------------------------------------------------------------------
        % erstelle Netz
%         [AC,AR] = meshgrid(ac,ar);
        [AC,AR] = ndgrid(ac,ar);
        
        % -------------------------------------------------------------------------
        % Interpolation von zwischenwerten
        
        YIA = 0;   % Eigentlich nicht Null aber Wert wird nicht gebraucht
        
        YIIA = griddedInterpolant(AC, AR, f2a', obj.InterpMethod);        
        
        YIIIA = 0;
        
        YIB = 0;   % Eigentlich nicht Null aber Wert wird nicht gebraucht
        
        YIIB = 0;   % Eigentlich nicht Null aber Wert wird nicht gebraucht
        
        YIIIB = griddedInterpolant(AC, AR, f3b', obj.InterpMethod);
        
    end % Ende Funktion
    
end % Ende Methoden

% -------------------------------------------------------------------------
% STATISCHE METHODEN
% ------------------------------------------------------------------------- 
methods (Static)    
    % ... bestimme schließspannung für einfaches ramberg osgood
    function scl = closesig(s,e,ecl,E,Ks,ns)
        % -----------------------------------------------------------------
        % Berechnung Rissschließdehnung aus RO Parameter
        % !!! mit verdopplungsgesetz auf halbästen
        % INPUT:
        %      s       - Spannung
        %      e       - Dehnung
        %      ecl     - Schließdehnung
        %      E       - E-Modul
        %    Ks,ns     - Parameter Ram. Osgood
        % OUTPUT:
        %     scl      - Schließspannung
        % -----------------------------------------------------------------
        
        % Sartwert
        scl = s - E*(e-ecl);
        % Fehler
        f = e - (s-scl)/E - 2 *((s-scl)/(2*Ks))^(1/ns)- ecl;
        % Newton
        while abs(f) > 1e-7
            % Ableitung
            df = 1/E - 2/ns * ((s-scl)/(2*Ks))^(1/ns-1)*(-1/(2*Ks));
            % neuer Wert
            scl = scl - 1/df * f;
            % Fehler
            f = e - (s-scl)/E - 2 *((s-scl)/(2*Ks))^(1/ns)- ecl;
        end
    end
       
    % ... vorzeichenbehaftete vergleichsspannung
    function [sv_vz, sw_vz] = mises_vorzeichen(DATA,imax)
        % Berechnet die Vorzeichenbehaftete von Mises Vergleichsspannungen
        % (vorzeichen nach der Normalspannungen auf der kritischen Ebene)
        % und die wirksame von Mises Vergleichsspannung in der kritischen Ebene
        %
        %
        % INPUT:
        %   DATA      -> Spannungs- und Dehnungswerte
        %   imax      -> Spaltenindex in der die maximale Normalspannung
        %                aufgetreteen ist
        %
        %
        % OUTPUT:
        % sv_max_vz   -> maximale Vorzeichenbehaftete von Mises Vergleichsspannung
        % sw_max_vz   -> maximale wirksame von mises Vergleichsspannung
        %
        %__________________________________________________________________________
        
        % Spannungen aus lesen
        sxx = DATA(1,imax);
        syy = DATA(2,imax);
        szz = DATA(3,imax);
        sxy = DATA(4,imax);
        syz = DATA(5,imax);
        sxz = DATA(6,imax);
        
        % von mises Vergleichsspannung und in der kritischen Ebene wirksame
        % Vergleichsspannung (mit vorzeichen aus normalspannung)
        if sxx >= 0
            sv_vz = sqrt( 0.5 * ((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + ...
                6 * ( sxy^2+syz^2+sxz^2 ) ) );
            sw_vz = sqrt( sxx^2 + 3 * (sxy^2 + sxz^2) );
        else
            sv_vz = -sqrt( 0.5 * ((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + ...
                6 * ( sxy^2+syz^2+sxz^2 ) ) );
            sw_vz = -sqrt( sxx^2 + 3 * (sxy^2 + sxz^2) );
        end
    end % Ende vorzeichenbehaftete Vergleichspannung
    
    % ... Rissöffnungsspannung
    function [sxop] = newman_max(sxo,sxu,sw,sv,sf,M_sig,ators)
        % Newman Verfahren zum bestimmen der Rissöffnungsspannung, für mehraxiale
        % Lasten angepasst von Döring
        %
        % INPUT:
        %       sxu    -> minimale Normalspannung in Ssp
        %       sxo    -> Smax Smax Suggar Smax
        %       sw    -> wirksame vorzeichenbehaftete Vergleichsspannung
        %       sv    -> vorzeichenbehaftete Vergleichsspannung
        %       sf    -> Fließspannung
        %      M_sig  -> Mittelspannungsempfindlichkeit zur Ermittlung des
        %                Parameters amitt
        %      astors -> Werkstoffparameter zur ermittlung der
        %                Schubspannungsempfindlichkeit
        %                (nach Vorschlag Hertl = 0 gesetzt)
        %
        % OUTPUT:
        %    sxop    -> Normalspannung bei der Riss geöffnet ist
        %
        %
        %__________________________________________________________________________
        
        
        if sxo > 0
            
            % Ermittle effektives R-Verhältniss
            % Normales R ( -oo < R <= 1 )
            R = sxu/sxo;
            
            % effektives R
            if R > 0
                Re = R;  % ( 0 < Re <= 1 )
            else
                Re = R * sv/sw; % ( 0 <= sv/sw <= 2 ); ( -oo < Rw <= 0 )
            end
            
            % Mittelspannungseinfluss (wie in FKM bei PJ)
            if M_sig < 0
                amitt = 0;
            else
                sxm = (sxo+sxu)/2;
                if sxm <= 1e-3 % eigentlich < 0, aber Toleranz
                    amitt = 0.4-M_sig/4;
                else
                    amitt = 0.47 -(1-1.5*M_sig)*(1+R)^(1+R+M_sig);
                end
            end
            
            
            % ermittle Konstanten
            A0 = 0.535 * cos(pi/2 * sw/sf) + amitt + ators * (1 - sv/sw);
            A1 = 0.344 * sw/sf + amitt;
            
            % Fallunterscheidung anhand effektiven R
            if Re > 0
                
                % weitere Konstanten
                A3 = 2*A0 + A1 - 1;
                A2 = 1 - A0 - A1 - A3;
                
                % Rissöffnungsspannung
                sxop = sxo * (A0 + A1*Re + A2*Re^2 + A3*Re^3);
                
            else
                
                % Rissöffnungsspannug
                sxop = sxo * (A0 + A1*Re);
                
            end
            
            % Abfangen sxop < sxu
            if sxop < sxu
                sxop = sxu;
            end
            
        else
            
            % Riss ist vollständig geschlossen
            sxop = sxo;
            
        end
    end % Ende Rissöffnungsspannung
    
    % ... Rissöffnungsdehnung
    function [exop,iop] = epsopen(DATA,iu,io,sxop,sxo)
        % Funktion zum auffinden der Rissöffnungsdehnung aus der
        % Rissöffnungsspannung, die sich bei einer einstufigen Belastung mit
        % konstanten amplituden einstellen würde.
        %
        % Annahmen:
        %   -> Suche auf Aufsteigendem Ast
        %   -> Bei mehreren Übereinstimmungen wird die kleinste Dehnung gewählt
        %   -> Zwischenwerte werden linear Interpoliert
        %
        %
        % INPUT:
        %    DATA   -> Spannungen und Dehnungen
        %    iu     -> index unten aufsteigender Ast
        %    io     -> index oben aufsteigender Ast
        %    sxop   -> Rissöffnungsspannung
        %    sxo    -> maximale Spannung im Ssp
        %
        %
        % OUTPUT:
        %   exop -> Rissöffnungsdehnung (linear interpolierter wert)
        %   iop  -> Index des ersten GEÖFFNETEN Punkts
        %
        % _________________________________________________________________________
        
        
        
        if sxop <= DATA(1,iu)
            % ... Riss ist immer offen
            iop = iu;
            exop = DATA(7,iop);
            
        elseif sxop > sxo
            % ... Riss ist immer zu
            iop = io;
            exop = DATA(7,iop);
            
        else
            % ... Riss teilweise offen
            
            % suche index indem öffnungsspannung überschritten wird
            iop = iu;
            while iop < io && DATA(1,iop) < sxop
                iop = iop + 1;
            end
            
            % lineare Interpolation der öffnungsdehnung
            if DATA(1,iop) - DATA(1,iop-1) ~= 0
                exop = (DATA(7,iop) - DATA(7,iop-1))/(DATA(1,iop) - DATA(1,iop-1)) ...
                    * (sxop - DATA(1,iop-1)) + DATA(7,iop-1);
            else
                exop = DATA(7,iop);
            end
        end
    end % Ende Rissöffnungsdehnung
    
    % ... bestimme schließspannung für beliebige hysteresen
    function [sxcl,icl] = sigclose(DATA,iu,io,excl,sxo,sxu)
        % Funktion zum auffinden der Rissschließspannugen aus der
        % Rissschließdehnung
        %
        % Annahmen:
        %   -> Suche auf Absteigendem Ast
        %   -> Zwischenwerte werden linear Interpoliert
        %
        %
        % INPUT:
        %    DATA   -> Spannungen und Dehnungen
        %    iu     -> index unten absteigender Ast (Endpunkt Halbast)
        %    io     -> index oben absteigender Ast (Startpunkt Halbast)
        %    excl   -> Rissschließdehnung
        %    sxo    -> maximale Spannung im Ssp
        %    sxu    -> minimale Spannung im Ssp
        %
        %
        % OUTPUT:
        %   sxcl -> Rissschließspannung (linear interpolierter wert)
        %   icl  -> Index des ersten GESCHLOSSENEN Punkts
        %
        % _________________________________________________________________________
        
        if excl >= DATA(7,io)
            % ... Riss immer zu
            icl = io;
            sxcl = sxo;
            
        elseif excl <= DATA(7,iu)
            % ... Riss immer auf
            icl = iu;
            sxcl = sxu;
            
        else
            % ... Riss teilweise offen
            
            % Suche index der Spannung zur Rissschließdehnung
            icl = io;
            while icl ~= iu && DATA(7,icl) > excl
                icl = icl + 1;
            end
            
            % ... Verhindere icl = 1
            %     if icl == 1
            %         icl = 2;
            %     end
            
            % lineare Interpolation
            sxcl = (DATA(1,icl)-DATA(1,icl-1))/(DATA(7,icl)-DATA(7,icl-1)) * ...
                (excl - DATA(7,icl-1)) + DATA(1,icl-1);
            
        end
    end % Ende Schließspannung
    
    % ... berechne Verzerrungsenergie
    function Weff = verzerrungsenergie(DATA,zs,ze,is,ie)
        % Funktion zum berechnen der Verzerrungsenergie der Hysteresehalbäste
        % Integration mit Trapezregel
        %
        %
        % INPUT:
        %  DATA   -> Spannungen und Dehnunge
        %  zs     -> Zeilenindex der Spannungen (zum unterscheiden der Moden)
        %  ze     -> Zeilenindes der Dehnungen (zum unterscheiden der Moden)
        %  is     -> Spaltenindex Start der Hysterese
        %  ie     -> Spaltenindex Ende der Hysterese (nicht am eigendlichen Ende
        %            sondern am schließzeitpunkt)
        %
        % OUTPUT:
        %  Weff   -> Verzerrungsenergie
        %
        %__________________________________________________________________________
        
        % Init
        Weff = 0;                    % Verzerrungsenergie
        sig0 = DATA(zs,is);          % Bezugsspannung
        
        % Trapezregel
        for n = is : ie-1
            
            Weff = Weff + 0.5 * (DATA(zs,n+1)-sig0+DATA(zs,n)-sig0) .* ...
                (DATA(ze,n+1)-DATA(ze,n));
            
        end % Ende Trapezregel Schleife
        
        % Weff = 0.5 * sum( ( DATA(zs,is+1:ie)+DATA(zs,is:ie-1)-2*sig0 ) .* ...
        %                   ( DATA(ze,is+1:ie)-DATA(ze,is:ie-1)) );
        
        % disp(Weff-Weff2);
    end % Ende Verzerrungsenergie
    
    % ... Schadensakkumulation 
    function [Dsum,Dakt] = DamAkk(P,k,ND,PD,Dsum)
        % Funktion zur Schadensakkumulation
        %
        % INPUT
        %  P     - Schädigungsparameter
        %  k     - Neigung P-WL
        % ND     - Dauerfestigkeit
        %  PD    - Dauerfestigkeit
        % Dsum   - Schädigung bisher
        %
        % OUTPUT:
        % Dsum   - aktualisierte Schadenssumme
        % Dakt   - aktueller Schädigungsbeitrag
        %__________________________________________________________________        
        if P > PD
            % ... Schädigung
            Dakt = 1/ND * (P/PD)^k;
            Dsum = Dsum + Dakt;
        else
            % ... keine Schädigung
            Dakt = 0;
        end
    end % Ende Schadensakkumulation
    
    % ... Simpson Integration
    function F = simpson(f,h,n)
        % Funktion zur Integration mit Simpson Regel
        %
        % INPUT:
        %  f       - Integranten an diskreten Punkten
        %  h       - Schrittweite
        %  n       - Anzahl Schritte
        %
        % OUTPUT
        %  F       - approx. Wert des Integrals F = int f(x) dx
        %__________________________________________________________________________
        
        % init mit start und endwerten
        F = f(1) + f(n);
        % schleife über alle anderen Werte
        for i = 2:n-1
            if mod(i,2) == 0 % gewicht = 2
                F = F + 2*f(i);
            else
                F = F +4*f(i);
            end
        end
        % alles Teilen
        F = h/3*F;
        
    end % Ende Simpson Integration

    % ... Interpolation zwischen Gestaltänderungsergie- und
    % Normalenhypothese
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
    end % Ende Interpolation zwischen Festigkeitshypothesen
    
end % Ende statische Methoden





end % Ende Klassendefinition