function [P] = read_RAINFLOW(fname)
% Funktion ließt Ergebnisse der Rainflowzählung aus .hcm Datei
%
% INPUT:
% fname     - Dateiname
%
% OUTPUT:
% P         - Schädigungsparameter
%            (1. Zeile = Mode nur bei Kurzriss)
%             1. Zeile = Schwingspiel
%             2. Zeile = Parameter
%             3. Zeile = Durchlaufzähler
%__________________________________________________________________________


% ... Anzahl Zeilen in Datei
[nl,nc] = GetNumberOfLine(fname);

% ... Speicher für OUTPUT
P = zeros(nc,nl-1);


% ... Öffne Datei mit Lesezugriff
fid = fopen(fname,'r');

% ... Lese Datei
status = 1;
z = 1;
while status
    % ... Lese nächste Zeile
    tline = fgetl(fid);
    if tline == -1
        % ... Lesen abbrechen
        status = 0;
    else
        % ... lese Werte
        [P(:,z),z] = readVal(tline,z,nc);
    end
end

% ... Schließe Datei
fclose(fid);

end % Ende Funktion

% ----------------------------------------------------------------------- %
% Hilfsfunktion
% ----------------------------------------------------------------------- %
function [p,z] = readVal(tline,z,nc)
    k = strfind(tline,';');
    p = zeros(nc,1);
    if length(k) == 2
        % ... SWT oder FS
        p(1,1) = str2double(tline(1:k(1)-1));
        p(2,1) = str2double(tline(k(1)+1:k(2)-1));
        p(3,1) = str2double(tline(k(2)+1:end));
        z = z + 1;
    elseif length(k) == 3
        % ... kurzriss
        p(1,1) = str2double(tline(1:k(1)-1));
        p(2,1) = str2double(tline(k(1)+1:k(2)-1));
        p(3,1) = str2double(tline(k(2)+1:k(3)-1));
        p(4,1) = str2double(tline(k(3)+1:end));
        z = z + 1;
    end

end

function [nl,nc] = GetNumberOfLine(fname)
    % ... Öffne Datei mit Lesezugriff
    fid = fopen(fname,'r');
    % ... Einfach Zeilen Zählen bis ende erreicht
    nl = 0;
    stat = 1;
    while stat
        tline = fgetl(fid);
        if tline == - 1
            stat = 0;
        else
            if nl > 0
                k = strfind(tline,';');
                nc = length(k)+1;
            end
            nl = nl + 1;
        end
    end
    % ... Schließe Datei
    fclose(fid);
end
