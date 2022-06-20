function [ZVAR0, EZVAR0] = ...
                          pseudostressapproach( ...
                          ZVAR0, EZVAR0, para, epara, esigpath, DEL, ...
                          ntens, ndi, material, numink, Outfile)
% Funktion führt Pseudo Stress Approach nach Kötten/Barkey/Socie aus
%
% AKTUELL NUR FÜR ESZ 
%
% INPUT:
%  ZVAR0, EZVAR0 -> Startwerte der zustandsvariabel Mat- und Str.Modell
%  para,epara    -> Parameter für mat und str modell
% esigpath       -> Verweiß auf Datei mit verlauf der lokalen pseudo
%                   Spannungen
% DE             -> Elastischer Nachgibigkeitstensor
% ntens, ndi     -> beschreiben Spannungszustand
% material       -> definiert welche materialfunktion aufgerufen wird (str)
% numink         -> Anzahl inkremente
%  M             -> Anzahl Backstresstensoren
% Outfile        -> Name der Datei in der die Ergebnisse geschrieben werden
%
% OUTPUT
%   Verlauf der Zustandsvariablen und ZV am Ende
%
% ------------------------------------------------------------------------
% Autor: Jan Kraft                                                        |
% Stand: April 2021                                                       |
% ------------------------------------------------------------------------

% Prüfe ob ebener Spannungszustand
if ntens ~= 3 && ndi ~= 2
    msg = 'Pseudo Stress Approach akt. nur für ESZ gedacht';
    error(msg)
end

% Öffne Inputdatei 
Infile = fopen(esigpath,'r');
% Startwerte
DATA0 = fread(Infile,[ntens+1,1],'double');
% Öffne Ergebnissdatei
fout = fopen(Outfile,'w');
% Buffer
BSize = 10000;

% materialfunktion
switch material
    case 'Chaboche'
        matfun = @chaboche;
    case 'OhnoWang'
        matfun = @ohnowang;
    case 'KarimOhno'
        matfun = @karimohno;
    case 'Jiang'
        matfun = @jiang;
    case 'Doring'
        matfun = @doring;
    case 'OWT'
        matfun = @OWT;
    otherwise
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end

% Ausführen Kerbsimulation (Hauptschleife über alle Werte in esigpath)
aktdat = 2;
tol = 1e-6;        % toleranz für iteration plastischer dehnungen
while aktdat <= numink
    space = min(BSize,numink-aktdat+1);                                    % Niemals mehr Speicher als Buffergröße freigeben
    % Einlesen Teilwerte aus pseudo Spannungsverlauf
    DATA = fread(Infile,[ntens+1,space],'double');                         % Einlesen Durchlaufzähler und pseudo Spannungen
    % Setzten Startwerte
    DATA = [DATA0,DATA];
    % Ermittle Pseudo Elastische Lösung
    ESIG = DATA(2:4,:);
    EEPS = DEL * ESIG;
    % Initialisieren der lokalen gesamt/plastischen Dehnungen
    EPSP = zeros(ntens,space+1);
    EPSP(:,1) = EZVAR0(ntens+1:2*ntens);
    
    % Initioalisieren Zustandsvariablen Materialmodell
    ZVAR = zeros(size(ZVAR0,1),space +1);
    ZVAR(:,1) = ZVAR0;
    
    % Initioalisieren Zustandsvariablen Strukturmodell
    EZVAR = zeros(size(EZVAR0,1),space+1);
    EZVAR(:,1) = EZVAR0;
    
    % Hauptschleife über alle Pseudo inkremente    
    for ii = 2 : space + 1
        
        % Pseudo Dehnungsinkrement
        dESIG = ESIG(:,ii) - ESIG(:,ii-1);
        
        % Integration Struckturmodell mit Pseudo Spannungsinkrement
        EZVAR0 = EZVAR(:,ii-1);
        EZVAR(:,ii) = matfun(ntens,ndi,dESIG,EZVAR0,0,epara);
        
        % Herauslesen der richtigen Gesamt-/Plastischendehnungen
        EPSP(:,ii) = EZVAR(ntens+1:2*ntens,ii);
        dEPSP0 = EPSP(:,ii) - EPSP(:,ii-1);
        
        % integration Materialmodell mit Spannungsinkrement
        fnorm = 1000;  % Abbruchbedingung für iteration
        iter = 0;      % zahler für while schleife
        dEPS = EEPS(:,ii)-EEPS(:,ii-1);    % erster versuch dehnungsinkrement
        ZVAR0 = ZVAR(:,ii-1);              % Zustandsvariablen zu beginn der Interation
        while fnorm > tol
            % Integration
            ZVAR(:,ii) = matfun(ntens,ndi,dEPS,ZVAR0,1,para);
            % berechne Fehler
            dEPSP = ZVAR(ntens+1:2*ntens,ii) - ZVAR0(ntens+1:2*ntens,:);
            f = dEPSP0 -dEPSP;
            % berechne Fehler norm
            fnorm = sqrt(f'*f);
            % berechne neues Dehungsinkrement falls keine konvergenz
            dEPS = dEPS + dEPSP0 - dEPSP;
            % inkrementieren iter
            iter = iter + 1;
            % verhindern von endlosschleifen
            if iter > 1000
                msg = ['In Lastschritt ' , num2str(ii), ' keine Konvergenz ',...
                    ' bei Bestimmung des plastischen Dehnungsinkrements'];
                error(msg)
            end
        end
        %     fprintf('Inkrement %i      Iterationen %i\n',ii,iter-1)
        % Ausführen des nicht ausgeführten Inkrements
        ZVAR(:,ii) = matfun(ntens,ndi,dEPS,ZVAR0,1,para);
        %     % Rücknahme des nicht ausgeführten Inkrements
        %     dEPS = dEPS - dEPSP0 + dEPSP;
        % inkrementieren der Dehnungen
        %     EPS(:,ii) = EPS(:,ii-1) + dEPS;
    end
    % pseudo Spannungen & DLZ am Ende
    DATA0 = DATA(:,space+1);
    % Zustandsvariablen am Ende
    ZVAR0 = ZVAR(:,space + 1);
    EZVAR0 = EZVAR(:,space + 1);
    % Herrauslesen der lokalen Größen je nach Material
    if aktdat == 2
        SIG = ZVAR(1:ntens,:);
        EPSP = ZVAR(ntens+1:2*ntens,:);
        EPS = DEL * SIG + EPSP;
        DLZ = DATA(1,aktdat-1:aktdat+space-1);
    else
        SIG = ZVAR(1:ntens,2:space+1);
        EPSP = ZVAR(ntens+1:2*ntens,2:space+1);
        EPS = DEL * SIG + EPSP;
        DLZ = DATA(1,2:space+1);
    end
    % Fehlende Dehnungskomponente
    [EPS, ~] = dehnungZZ(EPS,EPSP,para(2));
    % Rausschreiben Spannungen & Dehnungen 
    fwrite(fout,[DLZ;SIG;EPS],'double');
    % Inkrementiere Zeiger auf nächsten Index in ESIGALL
    aktdat = aktdat + space;                                               % Nächster Neuer Wert im nächsten Schleifendurchlauf
end % Ende Schleife über alle pseudo Spannungen

% Schließe Output & Input File
fclose(fout);
fclose(Infile);

end % Ende der Funktion