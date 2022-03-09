function [ZVAR0, EZVAR0] = devepsapproachV2( ...
                                   ZVAR0, EZVAR0, para, epara, esigpath, ...
                                   ntens, ndi, material, numink, M,eM, SIGO,...
                                   Outfile)
% Funktion führt deviatorischen Dehnungsapproach nach Hertl aus
%
% AKTUELL NUR FÜR ESZ 
%
% INPUT:
%  ZVAR0, EZVAR0 -> Startwerte der zustandsvariabel Mat- und Str.Modell
%                   (Schon umgestellt auf 3D)
%  para,epara    -> Parameter für mat und str modell
% esigpath       -> Verweiß auf Datei mit verlauf der lokalen pseudo
%                   Spannungen
% ntens, ndi     -> beschreiben Spannungszustand
% material       -> definiert welche materialfunktion aufgerufen wird (str)
% numink         -> Anzahl inkremente
%  M             -> Anzahl Backstresstensoren
%  eM            -> Anzahl Backstresstensoren Strukturmodell
% SIGO           -> Oberflächennormalspannung für hydrostatische Korrektur
%                   (0 im ESZ)
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
    msg = 'Pseudo Stress Approach akt. nur für ESZ gedacht';
    error(msg)
end

% Öffne Inputdatei 
Infile = fopen(esigpath,'r');
% Startwerte
DATA0 = fread(Infile,[ntens+1,1],'double');
% Öffne Ergebnissdatei
fout = fopen(Outfile,'w');
% fzvarout = fopen([Outfile(1:end-4),'.zvar'],'w');    % Zustandsvariablen (nur zum debuggen)
% Buffer
BSize = 10000;

% Elastizitätsparameter
E = para(1);
nu = para(2);
K = E/(3*(1-2*nu)); % kompressionsmodul
DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);

% ... Umstellen Zustandsvariablen auf 3D        
[ZVAR0] = umstellen_auf_3D(ZVAR0,material,M,nu,1);
[EZVAR0] = umstellen_auf_3D(EZVAR0,material,eM,nu,0);

% originale Tensoranzahl Speichern
ntensorig = ntens;
ndiorig = ndi;

% materialfunktion
switch material
    
    % Definiere Material Funktion
    case 'Chaboche'
        
        % definiere Materialfunktion
%         matfun = @chaboche; 
%         matfun = @chabocheRR2;
        matfun = @chabocheRR2_mex; 
        % Anzahl der Zustandsvariablen
        numzvar = (2+M)*6+2;
        enumzvar = (2+eM)*6+2;
        
        
    case 'OhnoWang'
        
        % Definiere Material Funktion
%         matfun = @ohnowang;
%         matfun = @ohnowangRR2;
        matfun = @ohnowangRR2_mex;
        % Anzahl der Zustandsvariablen
        numzvar = (2+M)*6+1;
        enumzvar = (2+eM)*6+1;
    
   case 'KarimOhno'
        
        % Definiere Material Funktion
%         matfun = @karimohnoRR2;
        matfun = @karimohnoRR2_mex;
        % Anzahl der Zustandsvariablen
        numzvar = (2+M)*6+1;
        enumzvar = (2+eM)*6+1;
        
    case 'Jiang'
        
        % Definiere Material Funktion
        matfun = @jiang;
        % Anzahl der Zustandsvariablen
        numzvar = (2+M)*6+2;
        enumzvar = (2+eM)*6+2;
        
    case 'Doring'
        % Definiere Material Funktion
        matfun = @doring;
        % Anzahl der Zustandsvariablen
        numzvar = (3+M+3)*7+1;
        enumzvar = (3+eM+3)*7+1;
        
    case 'OWT'
        % Definiere Material Funktion
        matfun = @OWT;
        % Anzahl der Zustandsvariablen
        numzvar = (3+M+(6+1)/2)*6+4;
        enumzvar = (3+eM+(6+1)/2)*6+4;
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
end

% andere tensoranzahl da in 3D leichter dev und hyd zu trennen geht
ntens = 6;
ndi = 3;

% Ausführen Kerbsimulation (Hauptschleife über alle Werte in esigpath)
aktdat = 2;
while aktdat <= numink
    space = min(BSize,numink-aktdat+1);                                    % Niemals mehr Speicher als Buffergröße freigeben
    % Einlesen Teilwerte aus pseudo Spannungsverlauf
    DATA = fread(Infile,[4,space],'double');                         % Einlesen Durchlaufzähler und pseudo Spannungen
    % Setzten Startwerte
    DATA = [DATA0,DATA];
    ESIG = DATA(2:4,:);
    
    % Startwert und Speicher der Zustandsvariablen
    ZVAR = zeros(numzvar,space+1);
    EZVAR = zeros(enumzvar,space+1);
    ZVAR(:,1) = ZVAR0;
    EZVAR(:,1) = EZVAR0;
    
    
    % Initialisieren der lokalen gesamt/plastischen Dehnungen und Spannungen
    % SIG = zeros(ntens,numink);
    EPS = zeros(ntens,space+1);
    EPS(:,1) = EZVAR(1:ntens,1);
    
    
    % Umstellen ESIG auf 3D ( nur ESZ )
    if ntensorig < 6
        ESIG = [ESIG(1,:);...
            ESIG(2,:);...
            zeros(1,space+1);...
            ESIG(3,:);...
            zeros(1,space+1);...
            zeros(1,space+1)];
    end
    
    
    % Hauptschleife über alle Pseudo inkremente
    for ii = 2 : space + 1
        
        % Pseudo Spannungsinkrement
        dESIG = ESIG(:,ii) - ESIG(:,ii-1);
        
        % Integration Struckturmodell mit Pseudo Spannungsinkrement
        EZVAR0 = EZVAR(:,ii-1);
        EZVAR(:,ii) = matfun(ntens,ndi,dESIG,EZVAR0,0,epara);
        
        % Herauslesen der Dehnungen/ des Dehungsinkrementes
        EPS(:,ii) = EZVAR(1:ntens,ii);
        
        % Inkrement der Dehnungen
        dEPS = EPS(:,ii) - EPS(:,ii-1);
        
        % Hydrostatisches (pseudo) Spannungsinkrement
        dESIGHYD = sum(dESIG(1:ndi))/3;
        
        % Nur Deviatorisches Dehnungsinkrement
        dEPSDEV = dEPS - dESIGHYD .* [1/(3*K);1/(3*K);1/(3*K);0;0;0];
        
        
        % Integration Materialmodell mit deviatorischem Dehungsinkrement
        % umstellen auf 3D, da hier leichter zu unterscheiden zwischen hyd und
        % dev
        ZVAR0 = ZVAR(:,ii-1);
        ZVAR(:,ii) = matfun(ntens,ndi,dEPSDEV,ZVAR0,1,para);
        
        
        % hydrostatisches Spannungsinkrement
        dSIGHYD = SIGO - ZVAR(3,ii);
        
        % Hydrostantische Korrektur
        ZVAR(1:ndi,ii) = ZVAR(1:ndi,ii) + dSIGHYD;
        
    end % Ende Schleife über Inkremente
    % pseudo Spannungen & DLZ am Ende
    DATA0 = DATA(:,space+1);
    % Zustandsvariablen am Ende
    ZVAR0 = ZVAR(:,space + 1);
    EZVAR0 = EZVAR(:,space + 1);
    % Herrauslesen der lokalen Größen je nach Material
    if aktdat == 2
        SIG = ZVAR([1 2 4],:);
        EPSP = ZVAR([7 8 10],:);
        EPS = DEL * SIG + EPSP;       
        DLZ = DATA(1,aktdat-1:aktdat+space-1);
    else
        SIG = ZVAR([1 2 4],2:space+1);
        EPSP = ZVAR([7 8 10],2:space+1);
        EPS = DEL * SIG + EPSP;
        DLZ = DATA(1,2:space+1);
    end
    % Fehlende Dehnungskomponente
    [EPS, ~] = dehnungZZ(EPS,EPSP,para(2));
    % Rausschreiben Spannungen & Dehnungen 
    fwrite(fout,[DLZ;SIG;EPS],'double');
    % Rausschreiben Zustandsvariablen (nur zum debuggen)
%     fwrite(fzvarout,[ZVAR(:,2:end);EZVAR(:,2:end)],'double');
    % Inkrementiere Zeiger auf nächsten Index in ESIGALL
    aktdat = aktdat + space;                                               % Nächster Neuer Wert im nächsten Schleifendurchlauf     
end % Ende Schleife über alle pseudo Spannungen

% Umstellen Zustandsvariablen auf ESZ
[ZVAR0] = umstellen_auf_ESZ(ZVAR0,material,M);
[EZVAR0] = umstellen_auf_ESZ(EZVAR0,material,eM);

% Schließe Output & Input File
fclose(fout);
fclose(Infile);
% fclose(fzvarout);

end