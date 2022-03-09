function superposition2(jobname,Lastpfad,c,ndl)
% Funktion erzeugen pseudo elastischer spannungsverlauf im Ebenen
% Spannungszustand und erzeugt Datei
%
% INPUT:
% jobname  - (str) Name der Rechnung
% Lastpfad - (str) Name der Lastdatei
% c        - Übertragungsfaktoren (double array, size: 3 x nkana)
% ndl      - Anzahl Durchläufe durch die Lastfolge, die in der
%           Kerbsimulation berechnet werden sollen
%
% OUTPUT:
% ESIG     - Verlauf pseudo elastische Spannungen 
%            (double array, size: 3 x ndata)
% DLZ      - Durchlaufzähler
% 
%
% -------------------------------------------------------------------------
% Autor: Jan Kraft (TU Darmstadt)
% Stand: April 21
% -------------------------------------------------------------------------

%% Init
BSize = 10000;
% ... abfangen fehler bei ndl
if ndl < 1, ndl = 1; end                                % weniger als 1 Durchlauf
if ~floor(ndl)==ndl, ndl = floor(ndl); end              % ndl kein integer

%% Bereite Lastfolgendatei & Outputdatei für pseudo Spannungen vor
% Lastdatei öffnen
fidLoad = fopen(Lastpfad,'r');
% Lese Header in Lastfolge:
tline = fgetl(fidLoad);
k = strfind(tline,':');
ndata = str2double(tline(k+1:end));
tline = fgetl(fidLoad);
k = strfind(tline,':');
nkana = str2double(tline(k+1:end));
tline = fgetl(fidLoad);
offset = ftell(fidLoad);                                                   % Read Position ausgehend von Start der Datei
% ESIG Datei öffnen
fidESIG = fopen(['00_Temp/',jobname,'.pth'],'w');
% Schreibe header 
header = ['    DLZ     ',...
          '            S11   ',...
          '          S22   ',...
          '          S12   \n'];
formspec = '%12.7f;%15.3f;%15.3f;%15.3f;\n';
fprintf(fidESIG,header);

%% Schreibe ESIG solange bis alle Durchläufe Erledigt sind
% Schleife über alle Durchläufe
dlz1 = 0;
for dl = 1 : ndl
    status = 1;                                                            % Abbruchbedingung einlesen
    % Lese Ganze Lastfolge
    counter = 1;                                                           % Zählt wie oft buffer gefüllt wird
    while status
        L = fscanf(fidLoad,'%f;%f;',[2,BSize]);                            % Lastfolgebuffer
        if ~isempty(L)
            % Berechne Anfahrt/Zwischenwerte 
            if counter == 1
                if dl == 1
                    L0 = L(:,1);                                           % Merke Startwerte
                    [LZW,nzw] = zwischen(L0);                              % Anfahrtweg
                else
                    [LZW,nzw] = zwischen(LE,L0);                           % Anfahrtweg
                end
            end
            % Berechne pseudo Elastische Spannungen
            ES = zeros(3,size(L,2)); 
            EZW = zeros(3,nzw); 
            for i = 1:nkana
                ES = ES + c(:,i) .* L(i,:);
                if nzw > 0
                    EZW = EZW + c(:,i) .* LZW(i,:);
                end
            end
            % Anhängen Zwischenwerte
            if counter == 1
                ES = [EZW,ES];
            end
            nsubdata = size(ES,2);
            % Berechne Durchlaufzähler
            dlz = dlz1 + 1/(nzw+ndata) : 1/(nzw+ndata) : dlz1 + nsubdata/(nzw+ndata);
            dlz1 = dlz1 + nsubdata/(nzw+ndata);
            % updates
            if size(L,2) < BSize
                status = 0;
            end
            counter = counter + 1;
            % Wegschreiben der Pseudo Elastischen Spannung
            fprintf(fidESIG,formspec,[dlz;ES]);
%             for i = 1:nsubdata
%                 fprintf(fidESIG,formspec,dlz(i),ES(:,i));
%             end
        else
            status = 0;
        end
    end % Ende Einlesen Lastfolge
    
    % Anhanägen Zwischenwerte
    LE = L(:,end);
    
    % Lastfolgedatei zurücksetzten
    fseek(fidLoad, offset, 'bof');    
    
end % Schleife über alle Durchläufe

%% Schließen aller offenen Dateien
fclose(fidLoad);
fclose(fidESIG);
end % Ende superposition




function [LZW,nzw] = zwischen(LE,L0)
% Berechnet Zwischenwerte zwischen zwei Durchläufen
% INPUT:
% LE    - Endwerte Lastfolge
% L0    - Startwerte der Lastfolge
% OUTPUT:
% LZW   - Zwischenwerte von L0 nach LE
% nzw   - Anzahl der Zwischenwerte
% -------------------------------------------------------------------------
if nargin == 1
    L0 = zeros(size(LE));
    t = 0:0.1:0.9;
else
    t = 0.1:0.1:0.9;
end
dL = LE-L0;
if sum(dL.*dL) <= 1e-3
    LZW = [];
    nzw = 0;
else
    LZW = t .* dL + L0;
    nzw = size(LZW,2);
end
end



