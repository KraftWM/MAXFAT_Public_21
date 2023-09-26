classdef Kerbsimulation < handle
% -------------------------------------------------------------------------
% Klassendefinition: Kerbsimulation 
% Als Zusammenfassung der Kerbsimulation. Speichert alle Optionen, steuert
% einzellne Berechnungsabläufe
%
%
%__________________________________________________________________________
% KONSTRUKTOR:
%    K = Kerbsimulation(jobname      - Name der Rechnung
%                       outpath      - Name Output Verzeichniss
%                       lastfolge,   - Matrix mit Lastfolge e R^(numkana x ndata)
%                       uefak,       - Matrix mit Werten der 
%                                      Übertragungsfaktoren
%                                      e R^(ntens x numkana)
%                       ndl,         - Anzahl Durchläufe durch die
%                                      Lastfolge  (int)
%                       Rm,          - Zugefestigkeit (double)
%                       E,           - Emodul
%                       nu,          - Querdehnzahl
%                       Kstrich      - Parameter Ramberg OSgood
%                       nstrich      - Parameter Ramberg Osgood
%                       varargin     - variabler Input
%                       'VarName1',var1,'VarName2',var2,...)
%
%  Es gibt Eingaben, die mindestens nötig sind um die Kerbsimulation
%  zu definieren, siehe oben. 
%  Für alle PUBLIC EIGENSCHAFTEN kann auch der Variablenname
%  und anschließend der Variablenwert übergeben werden. Die PUBLIC 
%  EIGENSCHAFTEN können auch nach Definition noch geändert werden. Das ist
%  über die Set-Methode der handle Klasse möglich
%    
%__________________________________________________________________________
%
% EIGENSCHAFTEN:
%  PUBLIC:
%
%   OPTIONEN:
%    jobname   -> Name der Rechnung (zum erzeugen von Dateien)
%    outpath   -> Name des Pfades zum erzeugen/lesen von Dateien
%    material  -> Name des materialmodells (string)
%                  "Chaboche"
%                  "KarimOhno" 
%                  "OhnoWang" (default)
%                  "Jiang"
%                  "OWT"
%   
% optpara,eoptpara -> Optionen zum Bestimmen der Material- & 
%                 Strukturparameter (wahl der stützstellen in plastischen 
%                                    Dehnungen)
%                   1 - äquidistant in Spannungen
%                   2 - geometrische reihe
%                   3 - verfahren simon moser (default)
%    
%   verfahren  -> Name des Näherungsverfahrens (string)
%
%                 Strukturfließflächen:
%
%                 "PseudoStress"
%                 "PseudoStress Lang" -> (default) Pseudo Stress aber schneller
%                 "PseudoStrain"
%                 "DevEps"
%
%                 Inkrementelle Verfahren:
%    
%                 "ESED"
%                 "Neuber"
%                 "ModNeuber"        -> Modifizierter Neuber
%                 "UniExp"
%                 
%
%  eindkerb    -> einachsiges Kerbnäherungsverfahren für Bauteilflieskurve
%                 "Neuber" =  default
%                 "ESED"
%                 "Seeger Beste"     -> Zusätzlicher Input Kp
%                 "Neuber Stern"     -> Zusätzlicher Input Kp
%                 "selbst"           -> Zusätzlicher Input BFK
%
%  nppara      -> Name des Nichtproportionalitätsparameters (string)
%                       "Riess"
%                       "Meggiolaro"
%                       "Gaier"
%                       "IQR"
%                       "Bishop"
%                       "Bolchoun"
%
%  NUMERISCHE FELDER:
% 
%  lastfolge   -> Lastfolge ( ein Durchlauf ) ) 
%   uefak      -> Übertragungsfaktoren
%   ndl        -> Anzahl zu berechnender Durchläufe durch die Lastfolge
%
%     Rm       -> Zugfestigkeit
%     M        -> Anzahl der zu verwendenden Backstresstensoren
%                default = 5
%    E         -> Emodul
%    nu        -> Querdehnzahl
%  Kstrich     -> Parameter zyklisches Ramberg Os. gesetz
%  nstrich     -> Parameter zyklisches Ramberg Os. gesetz
%   alpha      -> Materialparameter für Marquis/Socie ansatz
%    fk        -> Fließkurve, bietet die möglichkeit auf das Bauteil
%                 angepasste Fließkurven innerhalb des Plastizitätsmodells
%                 zu verwenden. Nicht zum
%                 speichern der mit 1d kerbnäherung ermittelten bfk
%                 1. Spalte = Spannungen   2. Spalte = plastische Dehnungen
%   bfk        -> Bauteilfließkurve, bietet die möglichkeit auf das Bauteil
%                 angepasste Fließkurven innerhalb der
%                 Struckturfließflächenansätze zu verwenden. Nicht zum
%                 speichern der mit 1d kerbnäherung ermittelten bfk
%                 1. Spalte = Spannungen   2. Spalte = plastische Dehnungen
%    Kp        -> Traglastformzahl ( falls Seeger Heuler oder Seeger Beste
%                 als 1d kerbnäherung für Bauteilfließkurve genommen wird) 
%    q         -> zum bestimmen der Materialparameter
%                 optoara = 1 -> q = maximale Spannung
%                 optpara = 2 -> q = Qutient geometrische Reihe
%                 optpara = 3 -> q = Anteil plastischer Dehnungen 
%                                            an gesamtdehnung bei
%                                            Fließbegin
%    eq        -> zum bestimmen der Strukturparameter, siehe q
%  ep_M        -> zum bestimmen der Materialparameter (nur bei optpara = 3)
%                 plastische dehunung nach der ideale plastizität gilt
%
%  eep_M       -> zum bestimmen der Strukturparameter(nur bei eoptpara = 3)
%                 plastische dehunung nach der ideale plastizität gilt
%   para       -> Parameter des Plastizitätsmodells (in PUBLIC damit sie 
%                 auch selbst gesetzt werden können)
%   epara      -> Parameter des Strukturmodells (in PUBLIC damit sie 
%                 auch selbst gesetzt werden können)
%   Cq         -> Parameter für Unified Expression
%   ZVAR       -> Zustandsvariablen Materialmodell
%  EZVAR       -> Zustandsvariablen Strukturmodell
%   REF        -> Referenzpunkte Inkrementelle Kerbnäherungen
%   chi        -> Ratchetting Parameter
%  echi        -> Pseudo Ratchetting Parameter
% npvhack      -> (bool) true  = NPV wird in SFF nur für Spannungen verwendet
%                                (default)
%                        false = NPV wird in SFF komplett verwendet
%
%  PROTECTED:
%   nkana      -> Anzahl Lastkanäle
%   ndata      -> Anzahl Datenpunkte
%    fnp       -> Bewertung des Lastpfades mit NP Kennzahl
%    Knp       -> Korrigierte Steifigkeit der Ramberg Osgood Gleichung nach
%                 Ansatz von Marquis und Socie
%
% maxSigMises  -> Maximale von Mises Vergleichsspannung in der Lastfolge
% esigpath     -> Verweis auf Datei in der Pseudoelastischer Spannungspfad
%                 gespeichert ist
%   ( ESIG      -> Pseudo Spannungsreihe Reihe)  -> wird jetzt aus Datei
%                                                         gelesen
%  (   DLZ       -> Durchlaufzähler
%  ( (ordnet Spannung einem Durchlauf/Zeitpunkt in der lastfolge zu) )
%  -> wird jetzt aus Datei gelesen
%                                                        
% 
%__________________________________________________________________________
% FUNKTIONEN:
% Normal:
%   bestimmeParameter -> Definiert Parameter für Kerbsimulation aus allen
%                        getroffenen Definitionen
%   kerbsimulation    -> Ausführen der Kerbsimulation mit gesetzten
%                        Definitionen
%   setVariableInput  -> Setzt variablen Input (protected)
%   NPKorrektur       -> Korrektur nach Marquis und Socie
% 
% Static:
%   checkLoadInput    -> Überprüft Input der Lastfolge
%    
%__________________________________________________________________________
% EXTRERNE FUNKTIONEN:
%   superposition   -> bestimmt pseudo elastische Lösung aus Lastfolge und
%                      Übertragungsfaktoren
%   ro2paraV2       -> bestimme Materialparameter aus ramberg osgood
%                      parametern
%   bfk2paraV2      -> bestimme Parameter aus Bauteilfließkurve
%   fk2para         -> unterfunktion aufgerufen in bfk2para & ro2para
%  init_matmodel    -> initialisieren startparameter der Materialmodelle
% neuberV4 ....     -> Inkrementelle Verfahren zur Kerbnäherung
% pseudo....        -> Strukturfließflächenansätze
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
      % Dateiverwaltung
      jobname                                                              % Name der Rechnung
      outpath                                                              % Name Verzeichniss zum speichern/laden von Dateien
      
      % Optionen
      material = 'OhnoWang';                                               % materialmodell
      verfahren = 'PseudoStress Lang';                                     % Name des verfahrens
      eindkerb = 'Neuber';                                                 % Einachsige Kerbnäherung für Bauteilfließkurve (Default = Neuber)
      nppara = 'None';                                                     % Welcher NP Parameter soll verwendet werden
      optpara = 3;                                                         % Methode zum bestimmen der Materialparameter
      eoptpara = NaN;                                                      % Methode zum bestimmen der Strukturparameter
      
      % Lastfolge
%       lastfolge {mustBeNumeric} = NaN;                                     % Lastfolge
      ndl {mustBeNumeric} = NaN;                                           % Anzahl Durchläufe durch die Lastfolge
      
      % Übertragungsfaktoren
      uefak {mustBeNumeric} = NaN;                                         % Übertragungsfaktoren
      
      % Zum Materialverhalten
      Rm {mustBeNumeric}                                                   % Zugfestigkeit
      E {mustBeNumeric}                                                    % Emodul
      nu {mustBeNumeric}                                                   % Querdehnzahl
      M {mustBeNumeric} = 10;                                              % Anzahl Backstresstensoren
	  eM {mustBeNumeric} = 10;                                             % Anzahl Backstresstensoren Strukturmodell
      Kstrich {mustBeNumeric}                                              % Para zyk Ram Os.
      nstrich {mustBeNumeric}                                              % Para zyk Ram Os.
      alpha {mustBeNumeric} = NaN;                                         % materialparameter für np hardening
      para {mustBeNumeric} = NaN;                                          % Parameter der funktion (definiert durch K', n' und Rm)
      epara {mustBeNumeric} = NaN;                                         % Parameter fürs struckturmodell
      Cq {mustBeNumeric} = NaN;                                            % Parameter für Unified Expression
      
      % Zum Bestimmen der Materialparameter
      q {mustBeNumeric} =  0.01;                                           % Parameter je nach Verfahren
      eq {mustBeNumeric} =  0.01;                                          % Parameter je nach Verfahren
      ep_M {mustBeNumeric} = 0.03;                                         % plastische Dehnung nach der ideale plasti herrscht
      eep_M {mustBeNumeric} = 0.03;                                        % plastische Dehnung nach der ideale plasti herrscht
      chi {mustBeNumeric} = NaN;                                           % Parameter für Ratchetting 
      echi {mustBeNumeric} = NaN;                                          % Pseudo Parameter für Ratchetting 
      
      % Zusätzliche Bauteilinformationen (Für Strukturfließflächenansätze)
      fk {mustBeNumeric} = NaN;                                            % Werkstofffließkurve
      bfk {mustBeNumeric} = NaN;                                           % Bauteilfließkurve
      Kp {mustBeNumeric} = NaN;                                            % Traglastformzahl
	  
	  % Zustandsvariablen und Referenzpunkte
	  ZVAR {mustBeNumeric} = NaN;                                          % Zustandsvariablen Materialmodell
	  EZVAR {mustBeNumeric} = NaN;                                         % Zustandsvariablen Strukturmodell
	  REF {mustBeNumeric} = NaN;                                           % Referenzpunkte inkrementelle Verfahren 

      % Berücksichtigen der Nichtproportionalen Verfestigung
      npvhack {mustBeNumericOrLogical} = true;

   end % Ende EIGENSCHAFTEN Public 
   
% PROTECTED
   properties (SetAccess = protected) % Größen werden automatisch gesetzt
      
      % Eigenschaften der Lastfolge 
      nkana {mustBeNumeric} = NaN;                                         % Anzahl Lastkanäle
      ndata {mustBeNumeric} = NaN;                                         % Anzahl Datenpunkte
      
      % materialmodell
      matfun                                                               % Funktion des materialmodells (gesetzt im Kontruktor)

      % nichtproportionalität
      fnp {mustBeNumeric}                                                  % Ergebnisse nichtproportionalitätsparameter
      Knp {mustBeNumeric}                                                  % Korrigierte Parameter RO Gleichung (nach Marquis und Socie)
      
      % Zeitreihen
      maxSigMises                                                          % Maximale von Mises Vergleichsspannung in der Lastfolge
	  maxSigModell                                                         % Spannungsgrenzwert gewähltes Materialmodell
	  maxESigModell														   % Spannungsgrenzwert gewähltes Strukturmodell
      esigpath                                                             % Verweis auf Pfad der Pseudoelastischen Lösung
%       ESIG  {mustBeNumeric}                                                % Pseudo Spannungsreihe
%       DLZ {mustBeNumeric}                                                  % Durchlaufzähler
      
            
   end % Ende EIGENSCHAFTEN protected 

% -------------------------------------------------------------------------
% KONSTRUKTOR
% ------------------------------------------------------------------------- 
   methods
       function obj = Kerbsimulation(jobname,outpath,...
                                     lastfolge,uefak,ndl,...
                                     Rm,E,nu,...
                                     Kstrich,nstrich,...
                                     varargin)
       %-------------------------------------------------------------------
       % Konstruktor der Klasse
       % INPUT:
       %    lastpfad        - Lastfolge e R^(numkana x numdata)
       %    uefak           - Übertragungsfaktoren e R^(3xnumkana)
       %    ndl             - Anzahl der Druchläufe
       %    Rm,E,nu         - statisches Materialverhalten
       %    Kstrich,nstrich - zyklisches Materialverhalten
       %    varargin        - cellarray
       %                     ('VarName1',var1,'VarName2',var2,...)
       % OUTPUT:
       %    obj - objekt der Klasse
       %-------------------------------------------------------------------
       
       % ------------------------------------------------------------------
       % Umgang mit benötigtem Input
       % ... Dateiverwaltung
       obj.jobname = jobname;
       obj.outpath = outpath;
       
       % ... felder zur Lastfolge
       [~, obj.nkana] = obj.checkLoadInput(lastfolge,uefak,ndl);
%        obj.lastfolge = lastfolge;
       obj.uefak = uefak;
       obj.ndl = ndl;
       
       % ... pseudo Spannungen durch Superpostion
       [ESIG,DLZ] = superposition(lastfolge,uefak,ndl);   
       obj.ndata = size(ESIG,2);
       % ... felder des Materialverhaltens
       obj.Rm = Rm;
       obj.E = E;
       obj.nu = nu;
       obj.Kstrich = Kstrich;
       obj.nstrich = nstrich;
       % ... maximale von Mises Vergleichsspannung
       obj.maxSigMises = max(vonMisesSpannung(ESIG));
       % ... Schreibe *.pth Datei
       obj.esigpath = [obj.outpath,'/',obj.jobname,'.pth'];
       fileID = fopen(obj.esigpath,'w');
       fwrite(fileID,[DLZ;ESIG],'double');
       fclose(fileID);
       
%        obj.esigpath = [obj.outpath,'/',obj.jobname,'.pth'];
%        write_SIGEPS(obj.jobname,obj.outpath,DLZ,ESIG)
       
       % ------------------------------------------------------------------
       % Umgang mit variablem Input
       obj.setVariableInput(varargin{:});     
       
       % ------------------------------------------------------------------
       % Nichtproportionalitätskorrektur
       obj.NPkorrektur(ESIG);
       
       % ------------------------------------------------------------------
       % Bestimme Strukturparameter
       obj.checkParaOpt();
       % ... Parameter sind nicht selbst gesetzt
       if isnan(obj.epara)
           obj.epara = obj.bestimmeParameter(obj.verfahren);
       end 
       % ... Grenzwert
       % obj.maxESigModell = spannungsgrenzwert(obj.material,obj.epara,obj.eM);


       % ------------------------------------------------------------------
       % Bestimme Materialparameter       
       % ... Parameter sind nicht selbst gesetzt
       if isnan(obj.para)
           obj.para = obj.bestimmeParameter("werkstoff");
       end 
       % ... Grenzwert
       % obj.maxSigModell = spannungsgrenzwert(obj.material,obj.para,obj.M);       
       
       
       % ------------------------------------------------------------------
       % Setze Startwerte
       % ... Startwerte nicht selbst gesetzt
       if isnan(obj.ZVAR)
           obj.ZVAR = init_matmodel(obj.material, 3,...
                        obj.para,obj.M);
       end
       
       % ... Startwerte nicht selbst gesetzt
       if isnan(obj.EZVAR) && ~any(isnan(obj.epara))
               obj.EZVAR = init_matmodel(obj.material, 3,...
                        obj.epara,obj.eM);
       end
       
       % ... Startwerte nicht selbst gesetzt
       if isnan(obj.REF)
           if strcmp(obj.verfahren,'Neuber')
               obj.REF = zeros(3,5);
           elseif strcmp(obj.verfahren,'ModNeuber')
               obj.REF = zeros(3,5);
           elseif strcmp(obj.verfahren,'UniExp')
               obj.REF = zeros(3,5);
           elseif strcmp(obj.verfahren,'ESED')
               obj.REF = zeros(3,3);
           end
       end
       
       end % Ende Konstrutor
       
   end % Ende KONSTRUKTOREN Public

% -------------------------------------------------------------------------
% METHODEN SET & GET
% ------------------------------------------------------------------------- 
    methods (Access = protected)
        % ... setze variablen Input
        function setVariableInput(obj,varargin)
            % Setzt alle variablen properties
            % 
            % INPUT:
            %    variabler Input
            %
            % OUTPUT
            %    setzt alle Properties die als Variabler Input übergeben
            %    werden können
            % -------------------------------------------------------------
            % ... Anzahl an variblem input
            nvarin = size(varargin,2);
            % ... checke ob nvarin gerade, falls nicht hinzufügen dummy wert
            if ~mod(nvarin,2), varargin{nvarin+1} = 0; nvarin = nvarin + 1; end
            % ... setze variablen
            for i = 1:2:nvarin-1
                varname = varargin{1,i};
                varvalue = varargin{1,i+1};
                % ... Unterscheide alle Eigenschaften
                switch varname
                    case 'material' % ... materialmodell
                        validvalue = {'Chaboche','KarimOhno','OhnoWang','Jiang','OWT'};
                        if any(strcmp(varvalue,validvalue))
                            obj.material = varvalue;
                        else
                            msg = ['Angegebenes Materialmodell nicht' ...
                                , ' erkannt. OhnoWang wird gesetzt'];
                            warning(msg);
                        end
                    case 'verfahren' % ... Verfahren
                        validvalue = {'Neuber','ESED','ModNeuber','UniExp',...
                            'PseudoStress','PseudoStress Lang',...
                            'DevEps','PseudoStrain'};
                        if any(strcmp(varvalue,validvalue))
                            obj.verfahren = varvalue;
                        else
                            msg = ['Angegebenes Verfahren nicht' ...
                                , ' erkannt. PseudoStress Lang wird gesetzt'];
                            warning(msg);
                        end
                    case 'optpara'  % ... Optionen Parameterbestimmung
                        if varvalue == 2
                            obj.optpara = varvalue;
                        elseif varvalue == 3
                            obj.optpara = varvalue;
                        elseif varvalue == 1
                            obj.optpara = varvalue;
                        else
                            obj.optpara = 3;
                            msg = ['Option zur Bestimmung der Parameter',...
                                ' nicht erkannt'];
                            warning(msg);
                        end
                    case 'eoptpara'  % ... Optionen Parameterbestimmung
                        if varvalue == 2
                            obj.eoptpara = varvalue;
                        elseif varvalue == 3
                            obj.eoptpara = varvalue;
                        elseif varvalue == 1
                            obj.eoptpara = varvalue;
                        else
                            obj.eoptpara = 3;
                            msg = ['Option zur Bestimmung der Parameter',...
                                ' nicht erkannt'];
                            warning(msg);
                        end
                    case 'nppara' % ... Nichtproportionalitätsparameter
                        validvalue = {'None','Riess','Meggiolaro','Gaier',...
                            'IQR','Bishop','Bolcoun'};
                        if any(strcmp(varvalue,validvalue))
                            obj.nppara = varvalue;
                        else
                            msg = ['Angegebenen NP Parameter nicht' ...
                                , ' erkannt. Rechnung ohne NP Parameter'];
                            warning(msg);
                        end
                    case 'eindkerb' % ... Einachsige Kerbnäherung für BFK
                        validvalue = {'Neuber','Neuber Stern','Seeger Beste',...
                            'ESED','selbst'};
                        if any(strcmp(varvalue,validvalue))
                            obj.eindkerb = varvalue;
                        else
                            msg = ['Angegebenen Verfahren zum Bestimmen der ',...
                                ' Bateilfließkurve nicht erkannt.',...
                                ' Rechnung mit einachsigem Neuber'];
                            warning(msg);
                        end
                    case 'M' % ... Anzahl Backstresstensoren
                        if isnumeric(varvalue) && varvalue > 0 && ~mod(varvalue,1)
                            obj.M = varvalue;
                        else
                            msg = ['Angegebene Anzahl an Backstresstensoren', ...
                                ' konnte nicht gesetzt werden'];
                            warning(msg);
                        end
                    case 'eM' % ... Anzahl Backstresstensoren
                        if isnumeric(varvalue) && varvalue > 0 && ~mod(varvalue,1)
                            obj.eM = varvalue;
                        else
                            msg = ['Angegebene Anzahl an Backstresstensoren', ...
                                ' konnte nicht gesetzt werden'];
                            warning(msg);
                        end
                    case 'alpha' % ... Nichtproportionalitätskennzahl Verfestigung
                        if isnumeric(varvalue) && varvalue > 0
                            obj.alpha = varvalue;
                        else
                            msg = ['Angegebenes alpha konnte nicht gesetzt',...
                                ' werden. Alpha wird automatisch bestimmt'];
                            warning(msg);
                        end
                    case 'para' % ... Materialparameter (kein Inputcheck)
                        obj.para = varvalue;
                    case 'epara' % ... Parameter Strukturmodell (kein Inputcheck)
                        obj.epara = varvalue;
                    case 'bfk' % ... Bauteilfließkurve (kein Inputcheck)
                        obj.bfk = varvalue;
                    case 'fk' % ... Fließkurve (kein Inputcheck)
                        obj.fk = varvalue;
                    case 'Cq' % ... Parameter Unified Expression
                        if varvalue < 0
                            obj.Cq = 0;
                        elseif varvalue > 1
                            obj.Cq = 1;
                        else
                            obj.Cq = varvalue;
                        end
                    case 'Kp' % ... Traglastformzahl
                        if isnumeric(varvalue) && varvalue > 1
                            obj.Kp = varvalue;
                        else
                            obj.Kp = 1.0001;
                            msg = ['Angegebenes Kp konnte nicht gesetzt',...
                                ' werden. Kp wird auf 1.001 gesetzt'];
                            warning(msg);
                        end
                    case 'q' % ... Parameter zu Bestimmung der Matpara
                        if isnumeric(varvalue) && varvalue > 0
                            obj.q = varvalue;
                        else
                            obj.q = 0.01;
                            msg = ['Angegebenes q konnte nicht gesetzt',...
                                ' werden. q wird auf 0.01 gesetzt'];
                            warning(msg);
                        end
                    case 'eq' % ... Parameter zu Bestimmung der Matpara
                        if isnumeric(varvalue) && varvalue > 0
                            obj.eq = varvalue;
                        else
                            obj.eq = 0.01;
                            msg = ['Angegebenes q konnte nicht gesetzt',...
                                ' werden. q wird auf 0.01 gesetzt'];
                            warning(msg);
                        end
                    case 'ep_M' % ... Parameter zu Bestimmung der Matpara
                        if isnumeric(varvalue) && varvalue > 0
                            obj.ep_M = varvalue;
                        else
                            obj.ep_M = 0.03;
                            msg = ['Angegebenes ep_M konnte nicht gesetzt',...
                                ' werden. ep_M wird auf 0.03 gesetzt'];
                            warning(msg);
                        end
                    case 'eep_M' % ... Parameter zu Bestimmung der Matpara
                        if isnumeric(varvalue) && varvalue > 0
                            obj.eep_M = varvalue;
                        else
                            obj.eep_M = 0.03;
                            msg = ['Angegebenes ep_M konnte nicht gesetzt',...
                                ' werden. ep_M wird auf 0.03 gesetzt'];
                            warning(msg);
                        end
                    case 'ZVAR' % ... Zustandsvariablen
                        obj.ZVAR = varvalue;
                    case 'EZVAR' % ... Zustandsvariablen Strukturmodell
                        obj.EZVAR = varvalue;
                    case 'REF' % ... Referenzpunkte
                        obj.REF = varvalue;
                    case 'chi' % ... Referenzpunkte
                        obj.chi = varvalue;
                    case 'echi' % ... Referenzpunkte
                        obj.echi = varvalue;
                    case 'npvhack'
                        obj.npvhack = varvalue;
                    otherwise  % ... Warnung bei Eigenschaft nicht erkannt
                        
                        msg = ['Die Eigenschaft ', varname, ' wurde nicht ',...
                            'als Eigenschaft der Klasse Kerbsimulation erkannt'];
                        warning(msg)
                end
            end % Ende Schleife Über den variablen Input
        end % Ende Setze variablen Input
         
        % ... checkt optionen für Parameter bestimmung
        function checkParaOpt(obj)
            % Prüft Input zum bestimmen der Strukturparameter
            % Falls nicht explizit ein Verfahren zum bestimmen der
            % Strukutrparameter angegeben wird, wird das gleiche Verfahren
            % wie für das Materialmodell verwendet.
            % 
            %  
            % -------------------------------------------------------------
            
            % ... Prüfe Strukturparameter
            if isnan(obj.eoptpara)
                obj.eoptpara = obj.optpara;
                obj.eep_M = obj.ep_M;
                obj.eq = obj.q;
                if obj.eoptpara == 1
                    obj.eq = 1.6 * obj.q;
                end
            end
        end % ende checken parameterbestimmung
    end % Ende Set & Get protected
    
    methods 
        
        % ... NP Korrektur
        function NPkorrektur(obj,ESIG)
            % Führt Korrektur der zyklischen Steifigkeit K' aus,
            % je nach nichtproportionalitäts
            %
            % Bestimmt fnp
            %
            % Setzt alle noch nicht gesetzten Felder 
            %
            % INPUT:
            % ESIG  - Lastpfad
            % 
            % OUTPUT
            %  Es werden nur Felder gesetzt
            % -------------------------------------------------------------
            % ... bestimme Kennwert
            switch obj.nppara
                case 'Riess'
                    obj.fnp = RiessFNP(ESIG');
                case 'Meggiolaro'
                    obj.fnp = MeggiolaroFNP(ESIG');
                case 'Gaier'
                    obj.fnp = GaierFNP(ESIG');
                case 'IQR'
                    obj.fnp = IQR_FNP(ESIG');
                case 'Bishop'
                    obj.fnp = BishopFNP(ESIG');
                case 'Bolchoun'
                    obj.fnp = BolchounFNP(ESIG');
                case 'None'
                    obj.fnp = 0;
            end % Ende NP Para
            % ... bestimme alpha aus Zugfestigkeit (falls nicht selbst gesetzt)
            % Diss Riess Seite 46 Gl. 4.4
            if isnan(obj.alpha)
                rp02 = obj.Kstrich * (0.002)^obj.nstrich;
                obj.alpha = ( (obj.Rm - rp02)/obj.Rm - 0.1 )/ 0.8;
            end
            % ... Korrektur
            obj.Knp = obj.Kstrich * ( 1 + obj.alpha * obj.fnp);
        end % Ende NP Korrektur
        
    end
% -------------------------------------------------------------------------
% METHODEN
% ------------------------------------------------------------------------- 
    methods
		
        % ... automatische Parameterbestimmung
        function parameter = bestimmeParameter(obj,typ)
            % -------------------------------------------------------------
            % Bestimme Material- oder Strukturparameter
            % INPUT:
            % typ  - "werkstoff" für Materialmodell 
            %        "Name des Verfahrens" für Strukturmodell
            % OUTPUT:
            % parameter - Parameter des Plastizitätsmodells
            % -------------------------------------------------------------
            
            % ... Dummy Wert
            parameter = NaN;
            % ... Alle StrukturFF Ansaätze
            validvalue = {'PseudoStress', 'PseudoStress Lang', ...
                         'DevEps','PseudoStrain'};
            % ... Organieren
%             Kwerkstoff = obj.Knp;
%             Ksturkt = obj.Knp;
            % ... Bestimme Materialparameter
            if strcmp(typ,'werkstoff')    
                if isnan(obj.fk)
                    if strcmp(obj.verfahren,'PseudoStrain')
                        if obj.npvhack
                            % NPV im Materialmodell ausstellen
                            parameter = ro2paraV2(typ, ...
                                        obj.E, obj.nu, obj.Kstrich, obj.nstrich,...
                                        obj.M, obj.material, ...
                                        obj.optpara,obj.q,obj.ep_M,...
                                        obj.eindkerb,obj.Kp,NaN,obj.chi);
                        else
                            % mit NP im Materialmodell
                            parameter = ro2paraV2(typ, ...
                                        obj.E, obj.nu, obj.Knp, obj.nstrich,...
                                        obj.M, obj.material, ...
                                        obj.optpara,obj.q,obj.ep_M,...
                                        obj.eindkerb,obj.Kp,NaN,obj.chi);
                        end
                    else
                        parameter = ro2paraV2(typ, ...
                                    obj.E, obj.nu, obj.Knp, obj.nstrich,...
                                    obj.M, obj.material, ...
                                    obj.optpara,obj.q,obj.ep_M,...
                                    obj.eindkerb,obj.Kp,obj.epara(end),obj.chi);
                    end
                else
                    if strcmp(obj.verfahren,'PseudoStrain')
                        parameter = bfk2paraV2(typ,obj.fk, ...
                                              obj.M, obj.material, ...
                                              obj.E, obj.nu,...
                                              obj.optpara,obj.q,obj.ep_M,NaN,obj.chi);
                    else
                        parameter = bfk2paraV2(typ,obj.fk, ...
                                              obj.M, obj.material, ...
                                              obj.E, obj.nu,...
                                              obj.optpara,obj.q,obj.ep_M,obj.epara(end),obj.chi);
                    end
                end
                % ... Sonderbehandlung OWT Modell
%                 if strcmp(obj.material,'OWT')
%                     [Qnpmax,gamma_np] = MarquisSocie2para(obj.Knp,obj.nstrich,obj.alpha);
%                     Qnpmax = 450;
%                     gamma_np = 100;
%                     parameter(7+3*obj.M) = Qnpmax;
%                     parameter(4+3*obj.M) = gamma_np;
%                 end
                % ... Bestimme Spannungsgrenzwert
                obj.maxSigModell = spannungsgrenzwert(obj.material,parameter,obj.M);
            % ... Bestimme Strukturparameter
            elseif any(strcmp(typ,validvalue))
                if strcmp(obj.eindkerb,'selbst') && ~any(any(isnan(obj.bfk)))
                    if strcmp(obj.verfahren,'PseudoStrain')
                        parameter = bfk2paraV2(obj.verfahren,obj.bfk, ...
                                              obj.eM, obj.material, ...
                                              obj.E, obj.nu,...
                                              obj.eoptpara,obj.eq,obj.eep_M,NaN,obj.echi);
                    else
                        parameter = bfk2paraV2(obj.verfahren,obj.bfk, ...
                                              obj.eM, obj.material, ...
                                              obj.E, obj.nu,...
                                              obj.eoptpara,obj.eq,obj.eep_M,obj.para(end),obj.echi);
                    end
                else
                    if strcmp(obj.verfahren,'PseudoStrain')
                        parameter = ro2paraV2(obj.verfahren, ...
                                obj.E, obj.nu, obj.Knp, obj.nstrich,...
                                obj.eM, obj.material, ...
                                obj.eoptpara,obj.eq,obj.eep_M,...
                                obj.eindkerb,obj.Kp,NaN,obj.echi);
                    else 
                        if obj.npvhack 
                            % NPV IM Strukturmodell ausstellen
                            parameter = ro2paraV2(obj.verfahren, ...
                                    obj.E, obj.nu, obj.Kstrich, obj.nstrich,...
                                    obj.eM, obj.material, ...
                                    obj.eoptpara,obj.eq,obj.eep_M,...
                                    obj.eindkerb,obj.Kp,obj.para(end),obj.echi);
                        else
                            % mit NPV im Strukturmodell
                            parameter = ro2paraV2(obj.verfahren, ...
                                    obj.E, obj.nu, obj.Knp, obj.nstrich,...
                                    obj.eM, obj.material, ...
                                    obj.eoptpara,obj.eq,obj.eep_M,...
                                    obj.eindkerb,obj.Kp,obj.para(end),obj.echi);
                        end
                    end
                end
                % ... Sonderbehandlung OWT Modell
%                 if strcmp(obj.material,'OWT')
%                     [Qnpmax,gamma_np] = MarquisSocie2para(obj.Knp,obj.nstrich,obj.alpha);
%                     Qnpmax = 450;
%                     gamma_np = 100;
%                     parameter(7+3*obj.M) = Qnpmax;
%                     parameter(4+3*obj.M) = gamma_np;
%                 end
                % ... Bestimme Spannungsgrenzwert
                obj.maxESigModell = spannungsgrenzwert(obj.material,parameter,obj.eM);
            end
        end % Ende bestimmeParameter
        
        % ... Kerbsimulation
        function [Outfile,ZVAR1,EZVAR1,REF1] = kerbsimulation(obj,...
                                                     varargin)
            % -------------------------------------------------------------
            % Ausführen in 3 Varianten:
            %       1. Nur Kerbsimulations Objekt wird übergeben
            %          -> Startwerte aller Größen werden selbst gesetzt
            %       2. Strukturfließfläche mit gegeben Startwerte
            %          -> varargin{1,1} = ZVAR0
            %          -> varargin{1,2} = EZVAR0
            %       3. Inkrementelle Verfahren mit gegebenen Startwerte
            %          -> varargin{1,1} = ZVAR0
            %          -> varargin{1,2} = REF0
            % INPUT:
            %   obj       - Objekt der Kerbsimulationsklasse
            %  varargin   - varibaler Input
            % OUTPUT:
            % Outfile - Verweis auf Datei mit Spannungs & Dehnungswerten
            % ZVAR1   - Zustandsvariablen am Ende der Kerbsimulation
            % EZVAR1  - Zustandsvariablen am Ende der Kerbsimulation
            %           (Strukturmodell)
            %  REF1   - Referenzwerte Umkehrpunkte am Ende der
            %           Kerbsimulation (inkrementelle Verfahren)
            % -------------------------------------------------------------
            % variabler Input:
            %  ZVAR0 EZVAR0   - Startwerte der Zustandsvariablen
            %  REF0           - Startwerte der Referenzpunkte für
            %                   Inkrementelle Verfahren
            % -------------------------------------------------------------  
            % -------------------------------------------------------------
            % ... Namer der Output Datei
            Outfile = [obj.outpath,'/',obj.jobname,'.sig'];
            
            % -------------------------------------------------------------
            % ... Unterscheide die Verschiedenen Varianten des Inputs
            strukff = {'PseudoStress', 'PseudoStress Lang', ...
                         'DevEps','PseudoStrain'};
            % 1. normal (Starte in 0)
            if nargin == 1 
                % ... bestimme Startwerte Materialmodell
                ZVAR0 = init_matmodel(obj.material, 3,...
                        obj.para,obj.M);
                % ... Strukturfließflächen
                if any(strcmp(obj.verfahren,strukff))
                    % ... Startwerte Strukturmodell
                    EZVAR0 = init_matmodel(obj.material, 3,...
                             obj.epara,obj.eM);
                    % ... Endwerte Referenzwerte 
                    REF1 = NaN;
                % ... Inkrementelle Verfahren
                else
                    % ... Endwerte Strukturmodell
                    EZVAR1 = NaN;
                    % ... Startwerte Referenzwerte
                    if strcmp(obj.verfahren,'Neuber')
                        REF0 = zeros(3,5);
                    elseif strcmp(obj.verfahren,'ModNeuber')
                        REF0 = zeros(3,5);
                    elseif strcmp(obj.verfahren,'UniExp')
                        REF0 = zeros(3,5);
                    elseif strcmp(obj.verfahren,'ESED')
                        REF0 = zeros(3,3);
                    end                    
                end
            % 2. Starte nicht in 0  
            % Ink Verfahren mit Startwert (!!! nicht getestet)
            elseif nargin == 3  && ~any(strcmp(obj.verfahren,strukff))
                ZVAR0 = varargin{1};
                REF0 = varargin{2};
            % Struk.FF Varfahren mit Startwert
            elseif nargin == 3  
                ZVAR0 = varargin{1};
                EZVAR0 = varargin{2};
            else
                msg = 'Falsche Anzahl Input';
                error(msg);
            end
            
            % -------------------------------------------------------------
            % ... Unterscheide Verfahren der Kerbsimulation
            switch obj.verfahren
%                 case 'PseudoStress'         % PseudoSpannungsansatz
%                     % ... Teste esigvmax < grenze Strukturmodell
%                     % ... berechne pseudo elastische Dehnung
%                     D = elast_nachgiebigkeit(obj.E,obj.nu,3,2);
%                     % ... aufrufen der Funktion 
%                     [ZVAR1,EZVAR1] = ...
%                      pseudostressapproach(...
%                      ZVAR0, EZVAR0, obj.para, obj.epara,...
%                      obj.esigpath, D,...
%                      3, 2, obj.material, ...
%                      obj.ndata, Outfile);
                case {'PseudoStress Lang' ,'PseudoStress'}  % PseudoSpannungsansatz schnell                    
                    % ... Teste ob vorgegebene Lastfolge Punkte außerhalb
                    % der Grenzen des Strukturmodells enthält
                    if obj.maxSigMises > obj.maxESigModell % spannungsgrenzwert(obj.material,obj.epara,obj.M)
                        % ... Setzte variablen input zum parameter anpassen
                        varinput = cell(1,3);
                        if ~strcmp(obj.eindkerb,'selbst') %   isnan(obj.bfk)
                            varinput{1} = obj.Knp;
                            varinput{2} = obj.nstrich;
                        else
                            varinput{1} = obj.bfk;
                        end
                        varinput{3} = obj.eq;
                        varinput{5} = obj.Kp;
                        % ... anpassen Strukturparameter
                        [obj.epara,obj.eq,obj.eep_M] = anpassenMaterialparameter(obj.maxSigMises,...
                                                              obj.verfahren,...
                                                              obj.eindkerb,...
                                                              obj.material,...
                                                              obj.M,...
                                                              obj.epara(1),obj.epara(2),...
                                                              obj.optpara,...
                                                              varinput{:});
                       % ... anpassen Materialparameter an
%                        obj.para = obj.bestimmeParameter('werkstoff');
                       % Strukturparameter
                    end
                    % ... Aufrufen des Pseudo Stress Approach nach lang
                    [ZVAR1, EZVAR1]  = ...
                     pseudostressapproachlangV2(ZVAR0,EZVAR0, obj.para, ...
                     obj.epara, obj.esigpath, ...
                     3, 2, obj.material, ...
                     obj.ndata, obj.M, obj.eM, Outfile);
                case 'PseudoStrain'         % PseudoDehnungsansatz
                    % ... berechne pseudo elastische Dehnung
                    D = elast_nachgiebigkeit(obj.E,obj.nu,3,2);
                    % ... Aufrufen des Pseudo Strain Approach
                    [ZVAR1,EZVAR1] = ...
                        pseudostrainapproach(ZVAR0, EZVAR0, obj.para, ...
                        obj.epara, obj.esigpath, D,...
                        3, 2, obj.material, ...
                        obj.ndata, Outfile);
                case 'DevEps'               % Deviatorischer Dehnungsansatz
                    % ... Aufrufen 
                    [ZVAR1,EZVAR1]  = ...
                        devepsapproachV2(ZVAR0, EZVAR0, obj.para, ...
                        obj.epara, obj.esigpath, ...
                        3, 2, obj.material, ...
                        obj.ndata, obj.M, obj.eM, 0, Outfile);
                case 'Neuber'               % max Neuber
                    % ... berechne pseudo elastische Dehnung
                    D = elast_nachgiebigkeit(obj.E,obj.nu,3,2);
                    % ... Aufrufen
                    [ZVAR1,REF1]  = ...
                    neuberV4(...
                    ZVAR0,REF0, obj.para, obj.esigpath,...
                    D, 3, 2, obj.material,...
                    obj.ndata,Outfile);
                case 'ESED'                 % max ESED
                    % ... berechne pseudo elastische Dehnung
                    D = elast_nachgiebigkeit(obj.E,obj.nu,3,2);
                    % ... aufrufen
                    [ZVAR1,REF1] = ...
                    esedV4(...
                    ZVAR0, REF0, obj.para, obj.esigpath,...
                    D, 3, 2, obj.material,obj.ndata,...
                    Outfile);
                case 'ModNeuber'            % Neuber mit Seeger Term
                    % ... berechne pseudo elastische Dehnung
                    D = elast_nachgiebigkeit(obj.E,obj.nu,3,2);
                    % ... aufrufen funktion
                    [ZVAR1,REF1] = ...
                    modneuberV4(...
                    ZVAR0, REF0, obj.para, obj.esigpath,...
                    D, 3, 2, obj.material,obj.ndata,...
                    Outfile);
                case 'UniExp'               % unified Expression
                    % ... berechne pseudo elastische Dehnung
                    D = elast_nachgiebigkeit(obj.E,obj.nu,3,2);
                    % ... berechne Dissipationskoef.
                    if isnan(obj.Cq)
                        obj.Cq = (1 - 2*obj.nstrich)/(1 - obj.nstrich);
                    end
                    % ... aufrufen funktion
                    [ZVAR1,REF1] = ...
                    uniexpV4(...
                    ZVAR0, REF0, obj.para, obj.Cq, obj.esigpath,...
                    D, 3, 2, obj.material,obj.ndata,...
                    Outfile);
            end                        
        end % Ende kerbsimulation
    end % Ende Methoden
    
% -------------------------------------------------------------------------
% STATISCHE METHODEN
% -------------------------------------------------------------------------    
    methods (Static)
        function [ndata, nkana] = checkLoadInput(L,c,ndl)
            % -------------------------------------------------------------
            % Bedingungen:
            %  - Lastfolge und Übertragungsfaktoren sind numerischer Input
            %  - gleiche Anzahl an Lastkanälen
            %  - 3 Übertragungsfaktoren für jeden Lastkanal (ESZ)
            %  - ndl ist positver integer
            % -------------------------------------------------------------
            
            
            % ... numerischer Input
            if ~isnumeric(L)
                eid = 'Input:notNumeric';
                msg = 'Die gegebene Lastfolge ist kein numerischer Input';
                throwAsCaller(MException(eid,msg))
            elseif ~isnumeric(c)
                eid = 'Input:notNumeric';
                msg = 'Die gegebenen Übertragungsfaktoren sind kein numerischer Input';
                throwAsCaller(MException(eid,msg))
            end
            
            % ... keine NaN Werte
            if sum(isnan(L),'all')
                eid = 'Input:isNaN';
                msg = 'Die gegebene Lastfolge enthält NaN Werte';
                throwAsCaller(MException(eid,msg))
            elseif sum(isnan(c),'all')
                eid = 'Input:isNaN';
                msg = 'Die gegebenen Übertragungsfakoren enthalten NaN Werte';
                throwAsCaller(MException(eid,msg))
            end
            
            % ... 3 Übertragungsfaktoren
            if size(c,1) ~= 3
                eid = 'Size:notEqual';
                msg = 'Falscher Spannungszustand in Übertragungsfaktoren';
                throwAsCaller(MException(eid,msg))
            end
            
            % ... Richtige Anzahl an Kanälen
            if size(c,2) ~= size(L,1)
                eid = 'Size:notEqual';
                msg = 'Anzahl Lastkanäle für Übertragungsfaktoren und Last passen nicht zusammen';
                throwAsCaller(MException(eid,msg))
            end
            
            % ... ndl ist positiver interger
            if ~isnumeric(ndl) || isnan(ndl) % kein numerischer Input
                eid = 'Input:notNumeric';
                msg = 'Die gegebene Anzahl an Durchläufen ist kein numerischer Input';
                throwAsCaller(MException(eid,msg))
            end
            
            % ... Alles richtig dann ermittle Arraygrößen
            nkana = size(L,1);  % Anzahl Kanäle
            ndata = size(L,2);  % Anzahl Datenpunkte
        end % Ende checkLoadInput
                
    end % Ende Statische Methoden
end % Ende Klassendefinition