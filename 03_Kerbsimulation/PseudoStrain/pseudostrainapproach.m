function [ZVAR0,EZVAR0] = ...
                          pseudostrainapproach( ...
                          ZVAR0, EZVAR0, para, epara, ...
                          esigpath, DE,...
                          ntens, ndi, material, numink, ...
                          Outfile)
% Funktion führt Pseudo Strain Approach nach Köttgen/Barkey/Socie aus
%
% AKTUELL NUR FÜR ESZ 
%
% INPUT:
%  ZVAR0, EZVAR0 -> Startwerte der zustandsvariabel Mat- und Str.Modell
%  para,epara    -> Parameter für mat und str modell
% esigpath       -> Verweiß auf Datei mit verlauf der lokalen pseudo
%                   Spannungen
% DE             -> Elastischer Nachgiebigkeitstensor
% ntens, ndi     -> beschreiben Spannungszustand
% material       -> definiert welche materialfunktion aufgerufen wird (str)
% numink         -> Anzahl der Inkremente
% Outfile        -> Name der Datei in der die Ergebnisse geschrieben werden
%
% OUTPUT
%   Zustandsvariablen am Ende
%
% ------------------------------------------------------------------------
% Autor: Jan Kraft                                                        |
% Stand: April 2021                                                       |
% ------------------------------------------------------------------------

% Prüfe ob ebener Spannungszustand
if ntens ~= 3 && ndi ~= 2
    msg = 'Pseudo Strain Approach akt. nur für ESZ gedacht';
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
%         matfun = @chaboche;
%         matfun = @chabocheRR2;
        matfun = @chabocheRR2_mex;
    case 'OhnoWang'
%         matfun = @ohnowang;
%         matfun = @ohnowangRR2;
        matfun = @ohnowangRR2_mex;
    case 'Jiang'
        matfun = @jiang;
    case 'KarimOhno'
%         matfun = @karimohno;
%         matfun = @karimohnoRR2;
        matfun = @karimohnoRR2_mex;
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
while aktdat <= numink
    space = min(BSize,numink-aktdat+1);                                    % Niemals mehr Speicher als Buffergröße freigeben
    % Einlesen Teilwerte aus pseudo Spannungsverlauf
    DATA = fread(Infile,[ntens+1,space],'double');                         % Einlesen Durchlaufzähler und pseudo Spannungen
    % Setzten Startwerte
    DATA = [DATA0,DATA];
    % Ermittle Pseudo Elastische Dehnung
    EEPS = DE * DATA(2:4,:);
    % Initialisieren der lokalen Spannungen
    SIG = zeros(ntens,space+1);
    SIG(:,1) = EZVAR0(1:ntens,1);
    
    % Initioalisieren Zustandsvariablen Materialmodell
    ZVAR = zeros(size(ZVAR0,1),space+1);
    ZVAR(:,1) = ZVAR0;
    
    % Initioalisieren Zustandsvariablen Strukturmodell
    EZVAR = zeros(size(EZVAR0,1),space+1);
    EZVAR(:,1) = EZVAR0;
    
    % Hauptschleife über alle Pseudo inkremente
    for ii = 2 : space + 1
        
        % Pseudo Dehnungsinkrement
        dEEPS = EEPS(:,ii) - EEPS(:,ii-1);
        
        % Integration Struckturmodell mit Pseudo Dehnungsinkrement
        EZVAR0 = EZVAR(:,ii-1);
        EZVAR(:,ii) = matfun(ntens,ndi,dEEPS,EZVAR0,1,epara);
        
        % Herauslesen der richtigen Spannungen
        SIG(:,ii) = EZVAR(1:ntens,ii);
        dSIG = SIG(:,ii) - SIG(:,ii-1);
        
        % integration Materialmodell mit Spannungsinkrement
        ZVAR0 = ZVAR(:,ii-1);
        ZVAR(:,ii) = matfun(ntens,ndi,dSIG,ZVAR0,0,para);
        
    end
    % pseudo Spannungen & DLZ am Ende
    DATA0 = DATA(:,space+1);
    % Zustandsvariablen am Ende
    ZVAR0 = ZVAR(:,space + 1);
    EZVAR0 = EZVAR(:,space + 1);
    % Herrauslesen der lokalen Größen je nach Material
    if aktdat == 2
        EPS = ZVAR(1:ntens,:);
        EPSP = ZVAR(ntens+1:2*ntens,:);
        DLZ = DATA(1,aktdat-1:aktdat+space-1);
    else
        SIG = SIG(:,2:space+1);
        EPS = ZVAR(1:ntens,2:space+1);
        EPSP = ZVAR(ntens+1:2*ntens,2:space+1);
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

end % Ende Funktion