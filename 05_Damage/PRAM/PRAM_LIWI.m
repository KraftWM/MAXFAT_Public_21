classdef PRAM_LIWI < handle
% -------------------------------------------------------------------------
% Klassendefinition: Schädigungsparameter PRAM
% Vergleichsspannung nach Alex Schmidt
%
%
%__________________________________________________________________________
% KONSTRUKTOR:
%    P = PRAM_LIWI(L,cf,ndl,...
%             E,Rm,...
%             Ks,ns,...
%             sf,ef,b,c,ND,...
%             Msig,nst,...
%             varargin);
%         L,cf,ndl       - Lastfolge, Übertragungsfaktoren und Anzahl
%                          Durchläufe
%         E,Rm           - statische Kennwerte
%         Ks,ns          - Ramber Osgood
%       sf,ef,b,c,ND     - Parameter Dehnungswöhlerlinie
%         Msig           - Mittelspannungsempfindlichkeit
%          nst           - Stützziffer
%           varargin     - variabler Input
%                       'VarName1',var1,'VarName2',var2,...)
% mit varargin können alle anderen zusätzlichen public Eigenschaften
% gesetzt werden
%__________________________________________________________________________
%
% EIGENSCHAFTEN:
% Public:
%  Ks,ns           - Ramber Osgood
% sf,ef,b,c        - Parameter Dehnungswöhlerlinie
%   ND             - Dauerfestigkeit Basis WL
%   E              - E Modul
%  Rm              - Zugfestigkeit
% nst              - Stützziffer
% Msig             - Mittelspannungsempfindlichkeit 
% miner            - option zum definieren der P-WL
%                    (0=elementar) (1=modifiziert) (sonst = orig)
%  ka              - Skalierung Normal- Schubbelastung in der
%                    Vergleichsspannungsdefinition
%                    (default = sqrt(3) -> Mises)
% method           - Amplitudenskalierung in der
%                    Vergleichsspannungsdefinition
%                     (default = 1 -> Preumont)
%  Kp              - Tranglastformzahl !!! Noch nicht berücksichtigt
%                     (default = 10000)
% Protected:
%   cc             - Counting Channel 
%                    !!! Eigentlich überflüssig, immer auf 1 lassen !!!
%  esigv           - Elastische Vergleichspannungen ( bestimmt aus
%                    Lastinformationen)
%  DLZ             - Fiktive Zeit
%   PN             - Logaritmierte Werte der P-WL
% Name             - Name des Parameters (nur zum generieren von Output)
%__________________________________________________________________________
% FUNKTIONEN:
% Normal:
%   pram        -> berechnet Schädigungsparameter
%   rainflow    -> rainflowzählung mit hcm
%   hcm         -> hcm zählung 
%   lebendauer  -> berechnet Lebensdauer
%   damage_akk  -> Schadensakkumulation
% 
% Static:
%   dwl2pwl          -> berechne PRAM - Wöhlerlinie aus Dehnungswöhlerlinie
% Vergleichsspannung -> Verlauf der Vergleichsspanung
%__________________________________________________________________________
% EXTRERNE FUNKTIONEN:
% Coin_liwi
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
      Msig {mustBeNumeric} = NaN;
      Rm {mustBeNumeric} = NaN;
      E {mustBeNumeric} = NaN;
      ND {mustBeNumeric} = NaN;
      miner {mustBeNumeric} = 3;
      Ks {mustBeNumeric} = NaN;
      ns {mustBeNumeric} = NaN;
      Kp {mustBeNumeric} = 10000;
      sf {mustBeNumeric} = NaN;
      ef {mustBeNumeric} = NaN;
      b {mustBeNumeric} = NaN;
      c {mustBeNumeric} = NaN;
      nst {mustBeNumeric} = 1;
      ka {mustBeNumeric} = sqrt(3);
      method {mustBeNumeric} = 1;
   end % Ende EIGENSCHAFTEN Public 
   
% PROTECTED
   properties (SetAccess = protected) % Größen werden automatisch gesetzt
      cc {mustBeNumeric} = 1;
      esigv {mustBeNumeric} = NaN;
      DLZ {mustBeNumeric} = NaN;
      PN {mustBeNumeric} = NaN;
      dauerfest = 1e20;         % definition für dauerfest
      Name = 'PRAM LIWI';            % Name des Parameters
   end % Ende EIGENSCHAFTEN protected 
   
% -------------------------------------------------------------------------
% KONSTRUKTOR
% ------------------------------------------------------------------------- 
methods
    function obj = PRAM_LIWI(L,cf,ndl,E,Rm,Ks,ns,sf,ef,b,c,ND,Msig,varargin)
        %------------------------------------------------------------------
        % Konstruktor der Klasse
        % INPUT:
        %         L,cf,ndl       - Lastfolge, Übertragungsfaktoren und 
        %                          Anzahl Durchläufe
        %      E                 - E Moduls
        %      Ks,ns             - Ramberg Osgood
        %    sf,ef,b,c,ND        - Parameter Dehnungswöhlerlinie
        %         Msig           - Mittelspannungsempfindlichkeit
        %          nst           - Stützziffer
        %           varargin     - variabler Input
        %                       'VarName1',var1,'VarName2',var2,...)
        % OUTPUT:
        %    obj - objekt der Klasse
        %------------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % Setze Felder
        obj.Rm = Rm;
        obj.Ks = Ks;
        obj.ns = ns;
        obj.sf = sf;
        obj.ef = ef;
        obj.b = b;
        obj.c = c;
        obj.E = E;
        obj.ND = ND;
        obj.Msig = Msig;
        
        % -----------------------------------------------------------------
        % Variabler Input
        if nargin >= 9
            nvarargin = length(varargin);
            for i = 2:2:nvarargin
                varname = varargin{i-1};
                varvalue = varargin{i};
                if strcmp(varname,'miner')
                    obj.miner = varvalue;
                elseif strcmp(varname,'ka')
                    obj.ka = varvalue;
                elseif strcmp(varname,'method')
                    obj.method = varvalue;
                elseif strcmp(varname,'nst')
                    obj.nst = varvalue;
                elseif strcmp(varname,'Kp')
                    obj.Kp = varvalue;
                end
            end
        end
       
        % -----------------------------------------------------------------
        % bestimme P-WL
        obj.PN = obj.dwl2pwl(E,sf,ef,b,c,ND,obj.nst,obj.miner);
        
        % -----------------------------------------------------------------
        % bestimme Verlauf Vergleichsspannung
        [obj.esigv,obj.DLZ] = obj.Vergleichsspannung(L,cf,ndl,obj.ka,obj.method);
        
       end % Ende Konstrutor
       
   end % Ende KONSTRUKTOREN Public   
   
% -------------------------------------------------------------------------
% METHODEN
% ------------------------------------------------------------------------- 
    methods
        % ... berechne Schädigungsparameter eines Schwingspiels
        function P = pram(obj,sna,snm,epsa)
            % Funktion zum berechnen des Schädigungsparameters nach 
            % Smith Watson Topper
            %
            % INPUT:
            %  DATA      -> Spannungen und Dehnungen und Druchläufe
            %   I        -> Index Start Hysterese
            %   J        -> Index UKP Hysterese
            %   K        -> Index Ende Hysterese
            %
            % OUTPUT:
            %  P         -> Schädigungsparameter des gezählten Ssp.
            %
            %__________________________________________________________________________
                        
            
            % ... Mittelspannungsempfindlichkeit
            if snm >= 0
                k = obj.Msig * (obj.Msig + 2);
            else
                k = obj.Msig/3 * (obj.Msig/3 + 2);
            end
            
            % -------------------------------------------------------------------------
            % Schädigungsparameter
            snmax = sna + k*snm;
            if snmax > 0
                P = sqrt(obj.E*snmax*epsa);
            else
                P = 0;
            end
        end % Ende Berechnung Schädigungsparameter
        
        % ... rainflowzählung
        function P = rainflow(obj)
            % -------------------------------------------------------------
            % Rainflowzählung mit Kerbnäherung
            %
            % OUTPUT:
            % P          - 1. Zeile Schwingspielzähler
            %              2. Zeile Schädigungsparameter
            %              3. Zeile Durchlaufzähler
            %              4. Zeile Spannungsamplitude
            %              5. Zeile Mittelspannung
            %              6. Zeile Dehnungsamplitude
            %              7. Zeile Mitteldehnung
            %
            % -------------------------------------------------------------
         
            % -------------------------------------------------------------------------
            % Bereitstellen Speicher
            Ndata = size(obj.esigv,1);                                                      % Anzahl Datenpunkte
            RES = zeros(3,Ndata);                                                      % reservierter Speicher für Residuum Mode I
            P = zeros(7,Ndata);                                                      % Schädigungsparameter (hier Amplituden da das nur ne testfunktion is)

            % -------------------------------------------------------------------------
            % Init Erstbelastungskurve & Hystereseast
            Sigma_IP=linspace(-5*obj.Rm,5*obj.Rm,150); 
            Sigma_IP = Sigma_IP';
            [Eps_EB,SigmaE_EB,~]=obj.zykSpannDehnNeuber(Sigma_IP,obj.E,obj.Ks,obj.ns);                 % Erstbelastung
            [Eps_H,SigmaE_H,~]=obj.HysteresenastNeuber(Sigma_IP,obj.E,obj.Ks,obj.ns);                   % Halbast (Schwingweiten)
            
            % -------------------------------------------------------------------------
            % Initialisierungen
            % ... Init Residuum (Zeiger, sig, eps)
            RES(1,1) = 1;                                                              % Startwert Residuum direkt auf 1
            RES(2,1) = interp1q(SigmaE_EB,Sigma_IP,obj.esigv(1));                             % Örtliche Spannungen
            RES(3,1) = interp1q(SigmaE_EB,Eps_EB,obj.esigv(1));                               % Örtliche Dehnungen
            % ... Init Zähler zum Speicher von Amplituden
            pz = 0;                                                                    % zeiger auf nächste frei Spalte in Pz
            % ... Zykluszähler
            counter = 0;                                 % Zykluszähler
            % ... Zeiger auf Residuum
            IR = 1;                                    % Zeiger auf Residuum
            IZ = 1;                                    % Zeiger auf Residuum
            
            
            % -------------------------------------------------------------------------
            % Rainflow Zählunf
            % ... Startindex
            K = 2;
            % ... Schleife über alle Werte
            while K <= Ndata

                % ... HCM aufrufen
                [IZ,IR,RES,counter,P,pz] = obj.hcm(K,IZ,IR,RES,counter,obj.esigv,P,pz,obj.DLZ,...
                                             Eps_EB,Sigma_IP,SigmaE_EB,...
                                             Eps_H,Sigma_IP,SigmaE_H);

                % ... Inkrementiere K
                K = K + 1;
            end
            
            % Freie Stellen in Output Leeren
            P = P(:,P(1,:) ~= 0);  % Aussortieren leere P Werte
                        
        end % Ende Rainflowzählung
        
        % ... hcm Algorithmus
        function [IZ,IR,RES,counter,P,pz] = hcm(obj,...
                K,IZ,IR,RES,counter,Data,P,pz,DLZ,...
                Eps_EB,Sigma_EB,SigmaE_EB,...
                Eps_H,Sigma_H,SigmaE_H)
            % Toleranzen
            tolM12 = 0.99;                             % Toleranz für das Erkennen vom Memory 1 und 2
            tolM3 = 1.01;                              % Toleranz für das Erkennen vom Memory 3
            sig = 0;                                   % Init dummy Wert
            eps = 0;                                   % Init dummy Wert
            
            % Abbruchbedingung
            weiter = 1;
            
            % ---------------------------------------------------------------------
            % Rainflow
            while weiter
                % 2
                if IZ > IR % Vergleich der Zeiger
                    
                    % ... letzte Werte aus Residuum lesen
                    I = RES(1,IZ-1);
                    J = RES(1,IZ);
                    
                    % ... Prüfe ob letzter Wert UKP ist
                    if (Data(K)-Data(J))*(Data(J)-Data(I)) >= 0
                        % ... kein UKP
                        IZ = IZ - 1;
                        weiter = 1;
                        %                 !! GOTO 2
                    else
                        % ... UKP
                        % ... Prüfe Schwingweite größer als die letzte
                        if abs(Data(K)-Data(J)) >= tolM12 * abs(Data(J)-Data(I))
                            % Örtlicher Zustand an UKP
                            sigI = RES(2,IZ-1);
                            epsI = RES(3,IZ-1);
                            sigJ = RES(2,IZ);
                            epsJ = RES(3,IZ);
                            % ... Schwingspiel gefunden
                            counter = counter + 1;
                            pz = pz + 1;
                            % ... Schwingspiele und Hysterese Schließindex
                            P(1,pz) = counter;
                            P(3,pz) = DLZ(K);
                            % ... Amplituden
                            siga = 0.5* abs(sigI-sigJ);
                            sigm = 0.5* (sigI+sigJ);
                            epsa = 0.5* abs(epsI-epsJ);
                            epsm = 0.5* (epsI+epsJ);
                            P(4,pz) = siga;
                            P(5,pz) = sigm;
                            P(6,pz) = epsa;
                            P(7,pz) = epsm;
                            P(2,pz) = obj.pram(siga,sigm,epsa);
                            IZ = IZ - 2;
                            % ... Stapel leer ?
                            if IZ >= IR
                                % ... nein, Memory 2
                                sig = RES(2,IZ) + interp1q(SigmaE_H,Sigma_H,Data(K)-Data(RES(1,IZ)));
                                eps = RES(3,IZ) + interp1q(SigmaE_H,Eps_H,Data(K)-Data(RES(1,IZ)));
                                %                         !! GOTO 2
                                weiter = 1;
                            else
                                % ... ja, Memory 1
                                sig = interp1q(SigmaE_EB,Sigma_EB,Data(K));
                                eps = interp1q(SigmaE_EB,Eps_EB,Data(K));
                                % ... Inkrementiere Zeiger
                                IZ = IZ + 1;
                                weiter = 0;
                            end
                        else
                            % ... UKP
                            sig = RES(2,IZ) + interp1q(SigmaE_H,Sigma_H,Data(K)-Data(J));
                            eps = RES(3,IZ) + interp1q(SigmaE_H,Eps_H,Data(K)-Data(J));
                            % ... kein Schwingspiel
                            weiter = 0;
                            IZ = IZ + 1;
                        end % Ende Verzweigung überprüfung der Schwingweiten
                    end % Ende Verzweigung UKP
                    
                else
                    % ... IZ <= IR (wird hier zsm behandelt)
                    % ... Einlesen Wert aus Residuum
                    J = RES(1,IZ);
                    weiter = 0;
                    % ... Prüfe UKP
                    if (Data(K)-Data(J))*Data(J) < 0
                        % ... Prüfe Memory 3
                        if abs(Data(K)) > tolM3 * abs(Data(J))
                            % ... UKP
                            sig = interp1q(SigmaE_EB,Sigma_EB,Data(K));
                            eps = interp1q(SigmaE_EB,Eps_EB,Data(K));
                            % ... Memory 3
                            IR = IR + 1;
                        else
                            % ... UKP
                            sig = RES(2,IZ) + interp1q(SigmaE_H,Sigma_H,Data(K)-Data(J));
                            eps = RES(3,IZ) + interp1q(SigmaE_H,Eps_H,Data(K)-Data(J));
                        end
                        IZ = IZ + 1;
                    end
                end % Ende Verzweigung Zeigervergleich
            end % Ende while Schleife
            
            % Speichern Residuum
            RES(1,IZ) = K;
            RES(2,IZ) = sig;
            RES(3,IZ) = eps;
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
            nssp = size(P,2);          % Anzahl Werte in P (=Anzahl Schwingspiele)
            Dsum = 0;                  % Schadenssumme
            Dlast = 0;                 % Schädigung letzter Durchlauf
            ndl = ceil(max(P(3,:)));   % maximale Anzahl an Durchläufen
            idam = 0;                  % Zeiger auf den Wert an dem Dsum = 1
            ilast = 1;                 % Zeiger auf ersten Wert des letzten Durchlaufs
            for i = 1: nssp
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
            if nssp > 0
                cyclast = P(1,nssp) - P(1,ilast);
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
                    %     Durchläufe
                    x = (1-Dsum)/Dlast;
                    DL = P(3,nssp) + x;
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
            %              1.Zeile logN
            %              2.Zeile logP
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
    end % Ende Methoden

% -------------------------------------------------------------------------
% STATISCHE METHODEN
% ------------------------------------------------------------------------- 
    methods (Static)
        % ... berechne P-WL aus D-WL
        function PN = dwl2pwl(E,sf,ef,b,c,ND,nst,miner)
            % ... faktor zum verschieben der WL !!! Ausgeschaltet
            fak = 1;% nst^(1/b);
            
            % ... WL
            npt = 4;                % Anzahl Punkte auf WL
            PN = zeros(2,npt);      % Speicher (1.Zeile = Ssp, 2.Zeile = P)
            for i = 1:npt
                % ... Punkt N
                N = ND^((i-1)/(npt-1));
                % ... Punkt P
                P =  sqrt(sf^2 * (2*N*fak)^(2*b)+E*sf*ef*(2*N*fak)^(b+c));
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
        
        % ... berechne Vergleichsspanung
        function [esigv,DLZ] = Vergleichsspannung(L,c,ndl,ka,method)
            % Größe des Lastsignals
            ndata = size(L,2);
            % Elastische Spannungen
            ESIG = repmat(c * L,1,ndl);
            % Fiktive Zeit
            DLZ = 1/ndata : 1/ndata : ndl;
            % Vergleichsspnnungen
            esigv = Coin_liwi(DLZ,ESIG,ka,method);
            % Filern Umkehrpunktfolge (1. und letzter Punkt werden erstmal
            % behalten=
            dL1 = [1;diff(esigv(1:end-1))];
            dL2 = [-1;diff(esigv(2:end))];
            idx = dL1 .* dL2 < 0;
            idx(length(esigv)) = 1;
            esigv = esigv(idx);
            DLZ = DLZ(idx);
            % Prüfe 1. Punkt !! Annahme Lastfolge startet in 0
            if sign(esigv(2)) == sign(esigv(1))         % Nur Filtern wenn Gleiches VZ
                if abs(esigv(2)) > abs(esigv(1))        % Wenn 1. Punk näher am Ursprung liegt
                    esigv(1) = [];
                    DLZ(1) = [];
                end
            end
        end
        
        % ... Hystereseast NeuberStern
        function [EpsP,SigmaE,EpsE]=HysteresenastNeuberStern(SigmaPo,E,K,n,Kp)
        % Berechnen der zyklische Spannungs-Dehnungs-Kurve
        % Als Input kommt die plastische Spannung und zyklische Parameter
        % daraus werden die zuhegörigen elastischen und elastisch-plastischen Dehnungen berechnet
        %
        % SigmaP:  elastisch-plastische Spannung
        % E: E-Modul
        % K: zyklischer Verfestigungskoeffizient
        % n: zyklischer Verfestigungsexponent
        %
        % Output:
        % EpsP: elastisch-plastisch Dehnung
        % EpsE: elastische Dehnung
        % SigmaE:  elastisch Spannung
        % -----------------------------------------------------------------
        
            SigmaP=abs(SigmaPo);
            
            eP=SigmaP./E+2.*(SigmaP./(2.*K)).^(1./n);
            N_P=SigmaP.*eP;
            
            
            SigmaEs=linspace(0,max(SigmaP)*5,1501);
            
            
            eS=(SigmaEs./Kp)./E+2.*((SigmaEs./Kp)./(2.*K)).^(1./n);
            N_E=SigmaEs.*Kp.*eS;
            
            SigmaE=sign(SigmaPo).*interp1q(N_E,SigmaEs,N_P);
            
            % zugehörige elastische Dehnung
            EpsE=SigmaE./E;
            
            
            % elastisch-plastische Dehnung
            EpsP=SigmaPo./E+sign(SigmaPo./E).*2.*(abs(SigmaPo)./(2.*K)).^(1./n);
            
        end % Ende Hystereseast Stern
        % ... Erstbelastung NeuberStern
        function [EpsP,SigmaE,EpsE]=zykSpannDehnNeuberStern(SigmaPo,E,K,n,Kp)
        % Berechnen der zyklische Spannungs-Dehnungs-Kurve
        % Als Input kommt die plastische Spannung und zyklische Parameter
        % daraus werden die zuhegörigen elastischen und elastisch-plastischen Dehnungen berechnet
        %
        % SigmaE:  elastisch-plastische Spannung
        % E: E-Modul
        % K: zyklischer Verfestigungskoeffizient
        % n: zyklischer Verfestigungsexponent
        %
        % Output:
        % EpsP: elastisch-plastisch Dehnung
        % EpsE: elastische Dehnung
        % SigmaE:  elastisch Spannung
            
            SigmaP=abs(SigmaPo);
            
            eP=SigmaP./E+(SigmaP./K).^(1./n);
            N_P=SigmaP.*eP;
            
            
            SigmaEs=linspace(0,max(SigmaP)*5,1501);
            
            
            eS=(SigmaEs./Kp)./E+((SigmaEs./Kp)./K).^(1./n);
            N_E=SigmaEs.*Kp.*eS;
            
            SigmaE=sign(SigmaPo).*interp1q(N_E,SigmaEs,N_P);
            
            % zugehörige elastische Dehnung
            EpsE=SigmaE./E;
            
            
            % elastisch-plastische Dehnung
            EpsP=SigmaPo./E+sign(SigmaPo./E).*(abs(SigmaPo)./K).^(1./n);
                        
        end
        % ... Hystereseast Neuber
        function [EpsP,SigmaE,EpsE]=HysteresenastNeuber(SigmaP,E,K,n)
        % Berechnen des Hysteresenasts
        % Als Input kommt die plastische Spannung und zyklische Parameter
        % daraus werden die zuhegörigen elastischen und elastisch-plastischen Dehnungen berechnet
        %
        % SigmaE:  elastisch-plastische Spannung
        % E: E-Modul
        % K: zyklischer Verfestigungskoeffizient
        % n: zyklischer Verfestigungsexponent
        %
        % Output:
        % EpsP: elastisch-plastisch Dehnung
        % EpsE: elastische Dehnung
        % SigmaE:  elastisch Spannung

            % elastische Spannung durch Neuberparabel aus plastischer Spannung
            SigmaE=2*sign(SigmaP).*(sqrt(abs(SigmaP/2).^2+E.*abs(SigmaP/2).*(abs(SigmaP/2)./K).^(1./n)));
        
            % zugehörige elastische Dehnung
            EpsE=SigmaE./E;
        
            % elastisch-plastische Dehnung
            EpsP=(SigmaP./E)+2*sign(SigmaP).*(abs(SigmaP)./(2.*K)).^(1./n);
        end
        % ... Erstbelastung Neuber
        function [EpsP,SigmaE,EpsE]=zykSpannDehnNeuber(SigmaP,E,K,n)
        % Berechnen der zyklische Spannungs-Dehnungs-Kurve
        % Als Input kommt die plastische Spannung und zyklische Parameter
        % daraus werden die zuhegörigen elastischen und elastisch-plastischen Dehnungen berechnet
        %
        % SigmaE:  elastisch-plastische Spannung
        % E: E-Modul
        % K: zyklischer Verfestigungskoeffizient
        % n: zyklischer Verfestigungsexponent
        %
        % Output:
        % EpsP: elastisch-plastisch Dehnung
        % EpsE: elastische Dehnung
        % SigmaE:  elastisch Spannung
            
            % elastische Spannung durch Neuberparabel aus plastischer Spannung
            SigmaE=sign(SigmaP).*sqrt(abs(SigmaP).^2+E.*abs(SigmaP).*(abs(SigmaP)./K).^(1./n));
            
            % zugehörige elastische Dehnung
            EpsE=SigmaE./E;
            
            % elastisch-plastische Dehnung
            EpsP=SigmaP./E+sign(SigmaP./E).*(abs(SigmaP)./K).^(1./n);
            
        end
    end % Ende statische Methoden

end % Ende Klassendefinition PRAM
