function [ZVAR0,REF] = modneuberV4( ...
                                   ZVAR0, REF0, para, esigpath, DEL, ...
                                   ntens, ndi, material, numink,...
                                   Outfile)
% Implentation einer Modifizierten Neuber Methode für Mehrachsige
% Kerbnäherungsverfahren.
%
% Idee, kombiniere Neuber mit Nebenbedingung für Umfangsspannung aus 
%       Hoffmannverfahren 
%
%  !!!!!!!  Implementierung für explizite Integration an normalen Punkten 
%           und implizite integration an umkehrpunkten             
%  !!!!!!!
%
% INPUT:
% ZVAR0     -> Startwerte der Zustandsvariablen, Zustandsvaraiblen für die
%              Materialmodelle wie bei dehnungsgesteuerter integration
% REF0      -> Startwerte der Referenzzustände
% para      -> Parameter Materialmodel
% esigpath  -> Verweiß auf Datei mit verlauf der lokalen pseudo
%              Spannungen
% DEL       -> Elastischer Nachgiebigkeitstensor
% ntens     -> Anzahl Tensorkomponenten
% ndi       -> Anzahl Hauptdiagonalelemente
% material  -> Definiert welches Materialmodell verwendet wird
% numink    -> Anzahl Inkremente der Lastfolge
% Outfile   -> Name der Datei in der die Ergebnisse geschrieben werden
%
% OUTPUT:
% ZVAR1  -> Zustandsvariablen am Ende
% REF1   -> Referenzzustände am Ende
%
% ------------------------------------------------------------------------
% Autor: Jan Kraft                                                        |
% Stand: April 2021                                                       |
% ------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Prüfe ob ebener Spannungszustand
if ntens ~= 3 && ndi ~= 2
    msg = 'Verfahren akt. nur für ESZ gedacht';
    error(msg)
end

% -------------------------------------------------------------------------
% Vorbereiten In-&Outputdateien
Infile = fopen(esigpath,'r');                                              % Öffne Inputdatei 
DATA0 = fread(Infile,[ntens+1,1],'double');                                % Startwerte
fout = fopen(Outfile,'w');                                                 % Öffne Ergebnissdatei
BSize = 10000;                                                             % Buffer

% -------------------------------------------------------------------------
% Definitionen für newton-verfahren
tol = 1e-16;                       % Abbruchkriterium
maxiter = 100;                      % Maximale Interation
alpha = 0.0;                       % Relaxationsparameter
% dummy = zeros(1,numink);            % Speichervaribale zum debuggen

% -------------------------------------------------------------------------
% materialfunktion
switch material
    case 'Chaboche'
        
        matfunexp = @chabocheenergyexp;
        matfunimp = @chaboche;
        
    case 'OhnoWang'
        
%         fprintf('noch ohnowang modell implementieren\n');
        matfunexp = @ohnowangenergy;
        matfunimp = @ohnowang;
        
    case 'KarimOhno'
        
%         fprintf('noch ohnowang modell implementieren\n');
        matfunexp = @karimohnoenergyexp;
        matfunimp = @karimohno;
        
    case 'Jiang'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end

% -------------------------------------------------------------------------
% Elastizitätskonstanten
E = para(1);                                                               % E-Modul
nu = para(2);                                                              % Querdehnzahl
Gmod = E / (2 * (1 + nu));
% CEL = elast_steifigkeit(E,nu,ntens,ndi);                                   % Steifigkeit
% DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);                                % Nachgiebigkeit
M = [2,-1,0;-1,2,0;0,0,3]./3;                                              % Abbildungen für Fließfläche
ML = [2,1,0;1,2,0;0,0,2];

% -------------------------------------------------------------------------
% Init Referenzzustände
% Annhame Lastfolge startet im Nullzustand -> Alle Refs mit 0
% inititalisieren
% Bei Neuber alle "realen" Spannungen an einem Referenzpunkt in Komponente 
% i speichern
% REF(3,5) = [REF_ESIG, REF_EEPS, REF_SIG, REF_EPSE, REF_EPSP]
REF = REF0;
dESIG0 = zeros(3,1);                                                       % Startwert Pseudospannungsinkrement

% Ausführen Kerbsimulation (Hauptschleife über alle Werte in esigpath)
aktdat = 2;
while aktdat <= numink
    % ---------------------------------------------------------------------
    % Pseudo elastische Werte einlesen
    space = min(BSize,numink-aktdat+1);                                    % Niemals mehr Speicher als Buffergröße freigeben
    % Einlesen Teilwerte aus pseudo Spannungsverlauf
    DATA = fread(Infile,[ntens+1,space],'double');                         % Einlesen Durchlaufzähler und pseudo Spannungen
    % Setzten Startwerte
    DATA = [DATA0,DATA];
    % Ermittle Pseudo Elastische Größem
    ESIG = DATA(2:4,:);
    EEPS = DEL * ESIG;
    
    % -------------------------------------------------------------------------
    % Finde ersten plastischen Zustand
    % ... Korrigiere erstes Inkrement im plastischen nur im ersten Durchlauf
    if sum(REF(:,1)) == 0 && aktdat == 2
        
        % Erster plastischer Zustand
        weiter = 1;
        jj = 1;
        r = para(end);
        while weiter && jj < space + 1
            s = M * ESIG(:,jj);
            F = s'*ML*s - 2/3*r^2;
            if F < 0
                jj = jj + 1;
            else
                jj = jj - 1;
                weiter = 0;
            end
        end
        
        % Automatisches Setzten von Zwischeninkrementen
        if jj == 1 % Erstes Inkrement plastisch -> Fehlermeldung
            %     msg = ['Fehler in Neuber. Erstes Inkrement plastisch, bitte mindestens',...
            %            ' ein elastisches Inkrement einfügen'];
            %     error(msg)
            msg = 'Erstes Inkrement ist plastisch, neues Inkrement wird berechnet';
            warning(msg)
            
            fak = 0.5;
            ESIG1 = fak * ESIG(:,2);
            s = M * ESIG1;
            F = s'*ML*s - 2/3*r^2;
            iter = 1;
            while F >= 0 && iter < 10
                fak = 0.5 * fak;
                ESIG1 = fak * ESIG(:,2);
                s = M * ESIG1;
                F = s'*ML*s - 2/3*r^2;
                iter = iter + 1;
            end
            ESIG = [ESIG(:,1),ESIG1,ESIG(:,2:end)];
            EEPS = [EEPS(:,1),DEL*ESIG1,EEPS(:,2:end)];
            DATA = [DATA(:,1),[0.5*(DATA(1,1)+DATA(1,2));ESIG1],DATA(:,2:end)];
            space = space + 1;
            jj = 2;
        end
    else
        jj = 1;
    end
    
    % -------------------------------------------------------------------------
    % Speicher für Zustandsvariablen
    ZVAR = zeros(size(ZVAR0,1),space + 1);                                 % Speicher
    
    % -------------------------------------------------------------------------
    % Startwerte der Zustandsvariablen
    if sum(REF(:,1)) == 0 && aktdat == 2
        % ... im 1. Durchlauf
        ZVAR(1:ntens,1:jj) = ESIG(:,1:jj);                                 % Speichern Elastische Spannungen
        ZVAR(ntens+1:end,1:jj) = ZVAR0(ntens+1:end).* ones(1,length(1:jj));% Speichern AB
    else
        % ... sonst
        ZVAR(:,1) = ZVAR0;
    end
    

    % -------------------------------------------------------------------------
    % Speicher Differenz aus soll Energieinkrement und ist Inkrement
    f = zeros(3,space + 1);
    
    % -------------------------------------------------------------------------
    % Hauptschleife über alle Inkremente der Lastfolge
    for ii = jj + 1 : space + 1
        
        % =====================================================================
        % Inkremente in Pseudo elastischen Größen
        dESIG = ESIG(:,ii) - ESIG(:,ii-1);
        dEEPS = EEPS(:,ii) - EEPS(:,ii-1);
        
        % =====================================================================
        % Prüfe auf Umkehrpunkt
        ukp = dESIG .* dESIG0 < 0;
        
        % =====================================================================
        % Prüfe auf erste Ecke
        eck = [false,false,false];
        for ee = 1 : 3
            if dESIG(ee) ~= 0 && dESIG0(ee) == 0 && ii > jj + 1
                eck(ee) = true;
            end
        end
        
        % =====================================================================
        % Speichern pseudo Spannungsinkrement falls nicht gleich 0
        idx = dESIG ~= 0;
        dESIG0(idx) = dESIG(idx);
        
        % =====================================================================
        % Update die Referenzzustände
        % REF(3,5) = [REF_ESIG, REF_EEPS, REF_SIG, REF_EPSE, REF_EPSP]
        EPSE = DEL * ZVAR(1:ntens,ii-1);                                       % Elastische Dehungen
        
        REF(ukp,1:4) = [ESIG(ukp,ii-1), EEPS(ukp,ii-1),ZVAR(ukp,ii-1),...
            EPSE(ukp)];
        for ee = 1 : ntens
            if ukp(ee)
                REF(ee,5) = ZVAR(ntens+ee,ii-1); % plastische Dehnungen
            end
        end
        
        % =====================================================================
        % Explitize Integration falls Inkrement kein Umkehrpunkt ist
        if sum(ukp) == 0 && sum(eck) == 0
            
            % Inkrement der totalen Verzerrungsenergiedichte des pseudo
            % Materials an ii - 1
            dETSED = (ESIG(:,ii-1) - REF(:,1)) .* dEEPS + ...
                (EEPS(:,ii-1) - REF(:,2)) .* dESIG;
            
            % PHAT
            EYY = EEPS(1,ii-1) - REF(1,2);
            EZZ = EEPS(2,ii-1) - REF(2,2);
            if EYY == 0
                EYY = 1e-8;
            end
            if EZZ == 0
                EZZ = (EEPS(2,ii-2)-REF(2,2))/(EEPS(1,ii-2)-REF(1,2))*1e-8;
                if isnan(EZZ)
                    EZZ = 0;
                end
            end
            PHAT = [ ZVAR(1,ii-1)-REF(1,3),     ZVAR(2,ii-1)-REF(2,3) , 0;...
                EZZ     ,                  -EYY     , 0;...
                0     ,                0          , ZVAR(3,ii-1)-REF(3,3)];
            
            % GHAT
            % GHAT = dESED/dSIG
            G11 = 2/E * ( ZVAR(1,ii-1)-REF(1,3) ) ...
                - 2 * nu / E * ( ZVAR(2,ii-1)-REF(2,3) ) ...
                + ( ZVAR(ntens+1,ii-1) - REF(1,5));
            G12 = 2/E * ( ZVAR(2,ii-1)-REF(2,3) ) ...
                - 2 * nu / E * ( ZVAR(1,ii-1)-REF(1,3) ) ...
                + ( ZVAR(ntens+2,ii-1) - REF(2,5));
            G21 = 1/E * ( EZZ + nu * EYY );
            G22 = - 1/E * ( nu * EZZ + EYY );
            G33 = 2/Gmod * ( ZVAR(3,ii-1)-REF(3,3) ) + ( ZVAR(ntens+3,ii-1) - REF(3,5));
            GHAT = [G11,G12,  0;...
                G21,G22,  0;...
                0,  0,G33];
            
            % Abfangen von Nullwerten in GHAT
            for ee = 1 : 3
                if GHAT(ee,ee) == 0
                    if dESIG(ee) == 0
                        GHAT(ee,ee) = 1;
                    else
                        GHAT(ee,ee) = dETSED(ee)/dESIG(ee);
                    end
                end
            end
            
            % Ableitung der Nebenbedingung
            dEPHI = [dETSED(1)+dETSED(2); ...
                0    ;...
                dETSED(3)];
            
            % explizite Integration des Materialmodells
            ZVAR(:,ii) = matfunexp(dEPHI, ZVAR(:,ii-1), para, GHAT, PHAT);
            
            % berechne tatsächliches Energieinkrement an ii-1
            dEPS = DEL * ZVAR(1:ntens,ii) +  ZVAR(ntens+1:2*ntens,ii) - ...
                DEL * ZVAR(1:ntens,ii-1) -  ZVAR(ntens+1:2*ntens,ii-1);
            SIG = ZVAR(1:ntens,ii-1);
            EPS = DEL * SIG + ZVAR(1+ntens:2*ntens,ii-1);
            dSIG = ZVAR(1:ntens,ii) - SIG;
            EPSREF = REF(:,4) + REF(:,5);
            dTSED = (SIG - REF(:,3)) .* dEPS + ...
                (EPS - EPSREF) .* dSIG;
            
            % berechne inkrement der nebenbedingungen an ii - 1
            phi2 = (EEPS(2,ii-1) - REF(2,2)) * dEPS(1) - ...
                (EEPS(1,ii-1) - REF(1,2)) * dEPS(2);
            dPHI = [dTSED(1) + dTSED(2); ...
                phi2; ...
                dTSED(3)];
            
            % =====================================================================
            % Impliziete Integration falls Inkrement ein Umkehrpunkt ist
        else
            
            % Inkrement der totalen Verzerrungsenergiedichte pseudo Material an
            % ii
            dETSED = (ESIG(:,ii) - REF(:,1)) .* dEEPS + ...
                (EEPS(:,ii) - REF(:,2)) .* dESIG;
            
            % Ableitung der Nebenbedingung
            dEPHI = [dETSED(1)+dETSED(2); ...
                0    ;...
                dETSED(3)];
            
            % Newton Raphson Verfahren
            [ZVAR(:,ii), dPHI] = newton_bisektion(...
                ZVAR(:,ii-1), dEPHI,...           % Zustand
                REF(:,3), REF(:,4) + REF(:,5),...% Refzustand
                DEL, para, matfunimp,...         % Material
                'ModNeuber',...                  % Definiert Näherungsverfahren
                maxiter,tol,alpha,...            % Iter optionen
                dEEPS,...                        % Startwert Dehungsinkrement
                ii, ...
                EEPS(:,ii), REF(:,2),dEEPS...
                );
            
            
        end % Ende Verzweigung Umkehrpunkte
        
        
        % =====================================================================
        % Abstandsquadrate soll - ist
        f(:,ii) = dEPHI - dPHI;
        
    end % Ende Schleife über Inkremente
    % pseudo Spannungen & DLZ am Ende
    DATA0 = DATA(:,space+1);
    % Zustandsvariablen am Ende
    ZVAR0 = ZVAR(:,space + 1);
    % Herrauslesen der lokalen Größen je nach Material
    if aktdat == 2
        SIG = ZVAR(1:ntens,:);
        EPSP = ZVAR(ntens+1:2*ntens,:);
        EPS = EPSP + DEL * SIG;
        DLZ = DATA(1,aktdat-1:aktdat+space-1);
    else
        SIG = ZVAR(1:ntens,2:space+1);
        EPSP = ZVAR(ntens+1:2*ntens,2:space+1);
        EPS = EPSP + DEL * SIG;
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