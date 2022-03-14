function [ZVAR1,EZVAR1] = ...
                          pseudostressapproachlangV2( ...
                          ZVAR0,EZVAR0, para, epara,esigpath, ...
                          ntens, ndi, material, numink, M,eM,...
                          Outfile)
% Funktion führt Pseudo Stress Approach nach Köttgen/Barkey/Socie mit dem
% Algo. von Lang et al. 
%
% AKTUELL NUR FÜR ESZ 
%
% INPUT:
%  ZVAR0, EZVAR0 -> Startwerte der Zustandsvariabel (schon richtig
%                   aufbereitet)
%  para,epara    -> Parameter für Mat.- und Str.-modell
% esigpath       -> Verweiß auf Datei mit verlauf der lokalen pseudo
%                   Spannungen
% ntens, ndi     -> beschreiben Spannungszustand
% material       -> definiert welche materialfunktion aufgerufen wird (str)
% numink         -> Anzahl inkremente
%  M             -> Anzahl Backstresstensoren
%  eM            -> Anzahl Backstresstensoren Strukturmodell
% Outfile        -> Name der Datei in der die Ergebnisse geschrieben werden
%
% OUTPUT
% Zustandsvariablen am Ende
%
% ------------------------------------------------------------------------
% Autor: Jan Kraft                                                        |
% Stand: April 2021                                                       |
% ------------------------------------------------------------------------
% =========================================================================
% Prüfe ob ebener Spannungszustand
if ntens ~= 3 && ndi ~= 2
    msg = 'Pseudo Stress Approach nach Lang akt. nur für ESZ gedacht';
    error(msg)
end
% =========================================================================
% Öffne Inputdatei & Lese Header
Infile = fopen(esigpath,'r');
% Startwerte
DATA0 = fread(Infile,[ntens+1,1],'double');

% =========================================================================
% Öffne Ergebnissdatein, schreibe Header & Definiere Outputformat
fout = fopen(Outfile,'w');
% fzvarout = fopen([Outfile(1:end-4),'.zvar'],'w');    % Zustandsvariablen (nur zum debuggen)

% =========================================================================
% Buffer
BSize = 10000;

% =========================================================================
% Zusammenfassen Zustandsvariablen für Material und Strukturmodell
ZVAR0 = pseudostresslang_ZV_zusammenfassen(...
                    ZVAR0,EZVAR0,material,M,eM);


% =========================================================================
% Definiere Materialfunktion und Zustandsvariablen
% Zuerst Zustandsvariablen des Struckturmodells (wie bei Spannungssteuerung)
% dann zusätzliche Variablen des Materialmodells 
switch material
    
    case 'Chaboche'
        
%         msg = 'Chabochemodell nicht implementiert für Ansatz nach Lang';
%         error(msg)
        % Speicher der Zustandsvariablen
%         nzvar = 7+6*M;      % Anzahl an Zustandsvariable, ohno Radien 
%       nzvar = (2+M+eM)*ntens + 1;
%         nzvar = 9+6*M;      % Anzahl an Zustandsvariable, mit Radien
        nzvar = (2+M+eM)*ntens + 3;
%         ZVAR = zeros(7+6*M,BSize);
        % Materialfunktion
%         matfun = @chabochelang;
%         matfun = @chabochelangRR2;
        matfun = @chabochelangRR2_mex;
%         matfun = @chabochelangRR2_2020b_mex;

    case 'OhnoWang'
        
        % Speicher der Zustandsvariablen
%         nzvar = 7+6*M;      % Anzahl an Zustandsvariable
        nzvar = (2+eM+M)*ntens+1;
%         ZVAR = zeros(7+6*M,BSize);
        % Materialfunktion
%         matfun = @ohnowanglang;
%         matfun = @ohnowanglangRR2;
        matfun = @ohnowanglangRR2_mex;
%         matfun = @ohnowanglangRR2_2020b_mex;

    case 'KarimOhno'
        
        % Speicher der Zustandsvariablen
%         nzvar = 7+6*M;      % Anzahl an Zustandsvariable
        nzvar = 1+ntens*(M+eM+2);      % Anzahl an Zustandsvariable
%         ZVAR = zeros(7+6*M,BSize);
        % Materialfunktion
%         matfun = @karimohnolang;
%         matfun = @karimohnolangRR2;
        matfun = @karimohnolangRR2_mex;
%         matfun = @karimohnolangRR2_2020b_mex;
    case 'Jiang'
        
        % Speicher der Zustandsvariablen
%         nzvar = 9+6*M;      % Anzahl an Zustandsvariable
        nzvar = (2+eM+M)*ntens+3;
%         ZVAR = zeros(9+6*M,BSize);
        % Materialfunktion
        matfun = @jianglang;
        
    case 'OWT'
        
        % Speicher der Zustandsvariablen
%         nzvar = (8+2*M)*3+7 ;      % Anzahl an Zustandsvariable
        nzvar = (5+M+eM+ntens)*ntens + 7;
        % Materialfunktion
        matfun = @OWTlang;    
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert für Ansatz nach Lang';
        error(msg)
end

% =========================================================================
% Initioalisieren Parameter für Materialmodell und struckturmodell
% para = [para, epara(3:end)];

% =========================================================================
% Ausführen Kerbsimulation (Hauptschleife über alle Werte in ESIGALL)
aktdat = 2;                                                                % Aktueller Zeilenindex in ESIGALL
while aktdat <= numink    
    space = min(BSize,numink-aktdat+1);                                    % Niemals mehr Speicher als Buffergröße freigeben
    % Einlesen Teilwerte aus pseudo Spannungsverlauf
    DATA = fread(Infile,[ntens+1,space],'double');                               % Einlesen Durchlaufzähler und pseudo Spannungen
    % Setzten Startwerte
    DATA = [DATA0,DATA];
    % Setze Startwerte der Zustandsvariablen
    ZVAR = zeros(nzvar,space+1);
    ZVAR(:,1) = ZVAR0;
    % Hauptschleife über alle Pseudo inkremente
    for ii = 2 : space+1

        % Pseudo Dehnungsinkrement
        dESIG = DATA(2:4,ii) - DATA(2:4,ii-1);

        % Integration Struckturmodell mit Pseudo Spannungsinkrement
        ZVAR0 = ZVAR(:,ii-1);
        ZVAR(:,ii) = matfun(dESIG,ZVAR0,para,epara);
%         ZVAR(:,ii) = matfun(dESIG,ZVAR0,para);

    end
    % pseudo Spannungen & DLZ am Ende
    DATA0 = DATA(:,space+1);
    % Zustandsvariablen am Ende
    ZVAR0 = ZVAR(:,space + 1);
    % Herrauslesen der lokalen Größen je nach Material
    if aktdat == 2
        [SIG, EPS, EPSP, ~, ~, ~, ZVAR1, EZVAR1] = ...
        pseudostresslang_ZV_auslesen(ZVAR,material,para,epara,M,eM,space+1);
        DLZ = DATA(1,aktdat-1:aktdat+space-1);
    else
        [SIG, EPS, EPSP, ~, ~, ~, ZVAR1, EZVAR1] = ...
        pseudostresslang_ZV_auslesen(ZVAR(:,2:space+1),material,para,epara,M,eM,space);
        DLZ = DATA(1,2:space+1);
    end
    % Fehlende Dehnungskomponente
    [EPS, ~] = dehnungZZ(EPS,EPSP,para(2));
    % Rausschreiben Spannungen & Dehnungen 
    fwrite(fout,[DLZ;SIG;EPS],'double');
    % Rausschreiben Zustandsvariablen (nur zum debuggen)
%     fwrite(fzvarout,ZVAR(:,2:end),'double');
    % Inkrementiere Zeiger auf nächsten Index in ESIGALL
    aktdat = aktdat + space;                                               % Nächster Neuer Wert im nächsten Schleifendurchlauf
end % Ende Schleife über alle pseudo Spannungen


% =========================================================================
% Schließe Output & Input File
fclose(fout);
fclose(Infile);
% fclose(fzvarout);
end