function [phic,psic,DLc,DL] = criticalplaneV3(...
                                      sigepsfile,ndata,...
                                      dphi,phimax,phimin,dpsi,psimax,psimin,...
                                      DMGs,...
                                      optdisplay,...
                                      optcritplane,...
                                      optrainflow,...
                                      optallhcm,...
                                      jobname,outpath)
%
% Funktion zum durchführen der kritischen Ebenen Schleife und
% Schädigungsrechnung
%
% Version benutzt Rainflowzählung aus Schädigungsparametern.
% Rainflowzählung mit Buffer
%
% -------------------------------------------------------------------------
% INPUT:
% Lokale Spannungen Dehnungen
%   ndata    -> (int) Anzahl Datenpunkte in sigepsfile (Anzahl Zeitpunkte
%               der Lstfolge)
% sigepsfile -> (str) Verweis auf Datei mit lokalem Lastpfad aus der
%               Kerbnäherung, enthält die folgenden Daten:
%   SIG      -> Lastfolge der Spannungen e R^(3,numink)
%   EPS      -> Lastfolge der Dehnungen e R^(3,numink)
%   EPSP     -> Lastfolge der plastischen Dehnungen e R^(3,numink)
%   DLZ      -> Ordnet Spannungen und Dehnungen den einzelnen Durchläufen
%               durch die Lastfolge zu
%
% Winkel Kritische Ebene
%   dphi     -> Inkrement des Winkels phi e R (Angabe in Grad)
%   dpsi     -> Inkrement des Winkels psi e R (Angabe in Grad)
% phimax,..  -> Maximal und Minimalwerte der Winkel für kritische Ebene
%              (Angabe in Grad)
%              (Default phi e [0, 180]; psi e [-90,90])
%
% Schädigungsmodell
% DMGs       -> cell array mit Objekten einer Schädigungsparameterklasse,
%                Der schädigungsparameter muss die Folgenden Funktionen
%                enthalten:
%                P = DMG.rainflow([sig;eps;DLZ]) -> Rainflow zum Identifizeiren und
%                                                   Berechnen der Schädigungsereignisse
%                DL = DMG.lebensdauer(P)         -> Schadensakkumulation 
%                (DMG.Name                       -> Name des Parameters
%                                          (nur für Display & Datei Output)
%
% Outputoptionen
% optdisplay     -> bool 1 = Display Ausgabe
%                        0 = keine Display Ausgabe
% optcritplane   -> bool 1 = Speichern der Ergebnisse Aller Ebenen in Datei
%                        0 = kein Speichern in Datei
% optrainflow    -> bool 1 = Speichern Rainflow Ergebniss der kritischen
%                            Ebene in Datei
%                        0 = kein Speichern in Datei
% opthcmall      -> bool 1 = Speichern Rainflow Ergbnisse ALLER Ebenen
%                        0 = nur speichern kritische Ebene
% jobname        -> str, name der Rechnung
% outpath        -> str, pfad für Output Dateien
%
% -------------------------------------------------------------------------
% OUTPUT:
% PHI     -> alle Winkel
% PSI     -> alle WInkel
% DL      -> Durchläufe 
% phic    -> phi der kritischen Ebene
% psic    -> psi der kritischen Ebene
% DLc     -> minimale Durchläufe
% _________________________________________________________________________
%
% ANMERKUNGEN
%
% Die Input Spannungen und Dehnungen befinden sich im Ebenen
% Spannungszustand. Es wird von lokalen Koordinatensystem {X,Y,Z} in der 
% Kerbe ausgegangen. Die Z-Achse steht dabei normal zur Bauteiloberfläche.
%
% Darstellung von Tensoren im globalen Koordinatensystem der Kerbe
%         sigXX              epsXX
%  sig =  sigYY       eps =  epsYY
%         sigXY             2epsXY
%
%
% Darstellung von Tensoren in gedrehten lokalen Koordinatensystemen (im
% Allgemeinen stellt sich 3D Spannungszustand ein). Spannungen und
% Dehnungen in den gedrehten Koordinatensystemen der kritischen Ebene
% {x,y,z}. Es wird dabei davon ausgegangen, dass die lokale x-Achse immer
% normal auf der betrachteten Schnittebene/dem Riss steht. die Spannung
% sigxx bezeichnet also die Normalspannung in der Schnittebene.
% -> y-Achse Bezeichnet Risslänge auf der Bauteiloberfläche
% -> z-Achse Bezeichnet Risslänge in Bauteiltiefe
% Bei dieser Wahl des lokalen Koordinatensystems lassen sich Normal- &
% Schubspannungen/-dehnungen direkt aus dem transformierten tensor ablesen
% es ist genau genommen keine weitere ebnen transformation mehr nötig.
%
%         sigxx              epsxx
%         sigyy              epsyy
%  sig =  sigzz       eps =  epszz
%         sigxy             2epsxy
%         sigyz             2epsyz
%         sigxz             2epsxz
%
% Definition der Winkel
% phi = Winkel zwischen Y,y => Rotation um Z e [0,pi]
% psi = Winkel zwischen Z,z => Rotation um y e [-pi/2,pi/2]
%
% _________________________________________________________________________

% ----------------------------------------------------------------------- %
% |                 Lade lokalen Lastpfad                               | %
% ----------------------------------------------------------------------- %
% [DLZ,SIG,EPS,~] = read_SIGEPS(sigepsfile);

% ----------------------------------------------------------------------- %
% |                 Umrechnene Winkelinkremente in rad                  | %
% ----------------------------------------------------------------------- %
dpsi = pi/180 * dpsi;
psimax = pi/180 * psimax;
psimin = pi/180 * psimin; 
dphi = pi/180 * dphi;
phimax = pi/180 * phimax;
phimin = pi/180 * phimin; 

% ----------------------------------------------------------------------- %
% |                 Speicher für Lebensdauerrechnung                    | %
% ----------------------------------------------------------------------- %
numdmg = length(DMGs);                                                     % Anzahl Schädigungsparameter
numpsi = ceil((psimax-psimin)/dpsi)+1;
numphi = ceil((phimax-phimin)/dphi)+1;
numwinkel = numpsi * numphi;
DL = zeros(numwinkel,2+numdmg);                                            % Speicher für output Variable

% ----------------------------------------------------------------------- %
% |           berechnen Dehnungskomponente in zz Richtung               | %
% ----------------------------------------------------------------------- %
% [EPS, ~] = dehnungZZ(EPS,EPSP,nu);


% ----------------------------------------------------------------------- %
% |                 Schleife über alle Ebenen                           | %
% ----------------------------------------------------------------------- %
% ... Init Variablen
zahler_planes = 1;   % Zähler variable der ebenen
DLc = 1e21 * ones(1,numdmg);        % minimale Durchläufe initial sau hoch setzten
Pcrit = struct();                  % Schädigungsparameter in kritischer Ebene
phic = zeros(1,numdmg);            % Winkel kritische Ebene
psic = zeros(1,numdmg);            % Winkel kritische Ebene

% ... Display Ausgabe
if optdisplay
    fprintf(' plane phi   psi   ');
    for i = 1 : numdmg
        sizename = length(DMGs{i}.Name);
        fprintf([' ',DMGs{i}.Name]);
        for j = 1:13-sizename
            fprintf(' ')
        end
    end
    fprintf('\n');
end
        
% ... Schleife
for phi = phimin : dphi : phimax % Drehung um Z

    for psi = psimin : dpsi : psimax % Drehung um y
        
        % ... Display Ausgabe
        if optdisplay
            fprintf('%6i%6.3f%6.3f',zahler_planes,phi,psi);
        end
        
        % ... Abspeichern Winkel
        DL(zahler_planes,1,:) = phi * 180/pi;
        DL(zahler_planes,2,:) = psi * 180/pi;
        
        % ... Transformation der Spannungen ins lokale Koordinatensystem
%         sig = transformstress(SIG,phi,psi);
        
        % ... Transformation der Dehnungen ins lokale Koordinatensystem
%         eps = transformstrain(EPS,phi,psi);
        
        % ... Transformation der plast. Dehnungen ins lokale Koordinatensystem
%         epsp = transformstrain(EPSP,phi,psi);
        
        % --------------------------------------------------------------- %
        % |                 Schädigungsrechnung                         | %
        % --------------------------------------------------------------- %
        for i = 1:numdmg
            % ... aktueller Schädigungsparameter
            DMG = DMGs{i};
            
            % ... Rainflowzählung und berechnen Schädigungsparameter
%             P = DMG.rainflow([sig;eps;DLZ]);
            P = DMG.rainflow(sigepsfile,ndata,phi,psi);
            if optallhcm && optrainflow
                write_RAINFLOW(jobname,DMG.Name,outpath,P,phi,psi);
            end
            
            % ... Lebendsdauer berechnen
            DL(zahler_planes,2+i) = DMG.lebensdauer(P);
            
            % ... merke kritische ebene
            if DL(zahler_planes,2+i) < DLc(i)
                DLc(i) = DL(zahler_planes,2+i);
                Pcrit.(DMG.Name) = P;
                phic(i) = phi;
                psic(i) = psi;
            end
            
            
            % ... Display Ausgabe
            if optdisplay
                fprintf('%14.6d',DL(zahler_planes,2+i));
            end
        end
        
        % ... Display Ausgabe
        if optdisplay
            fprintf('\n');
        end
        
        % ... Inkrementieren des ebenen zählers
        zahler_planes = zahler_planes + 1;
        
    end % Ende Schleife psi
    
end % Ende Schleife phi

% ... Dateiausgabe kritische Ebene
if optcritplane
    write_CRITPLANE(jobname,outpath,DMGs,DL);
end

% ... Dateiausgabe Rainflow kritische Ebene
if optrainflow && ~optallhcm
    for i = 1:numdmg
        write_RAINFLOW(jobname,DMGs{i}.Name,outpath,Pcrit.(DMGs{i}.Name),phic(i),psic(i))
    end
end

end % Ende Funktion