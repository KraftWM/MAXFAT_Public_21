function L = read_LASTFOLGE(fname)
% Funktion ließt Spannungs und Dehnungswerte aus .sig datei ein
%
% INPUT:
% fname     - Dateiname
%
% OUTPUT:
% L         - Lastfolge
%__________________________________________________________________________


% % ... Anzahl Zeilen in Datei
% [nl,nk] = GetNumberOfLine(fname);

% ... Öffne Datei mit Lesezugriff
fid = fopen(fname,'r');

% ... Anzahl Daten/ Zeilen
tline = fgetl(fid);
k = strfind(tline,':');
nl = str2double(tline(k+1:end));

% ... Anzahl Lastkanäle
tline = fgetl(fid);
k = strfind(tline,':');
nk = str2double(tline(k+1:end));


% ... Speicher für OUTPUT
L = zeros(nk,nl);

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
        [L(:,z),z] = readVal(tline,z,nk);
    end
end

% ... Schließe Datei
fclose(fid);

end % Ende Funktion


% ----------------------------------------------------------------------- %
% Hilfsfunktion
% ----------------------------------------------------------------------- %
function [l,z] = readVal(tline,z,nk)
    k = strfind(tline,';');
    l = zeros(nk,1);
    if ~isempty(k)
        l(1) = str2double(tline( 1:k(1)-1 ));
        for i = 2:length(k)
            l(i) = str2double(tline( k(i-1)+1:k(i)-1 ));
        end
        z = z + 1;
    end

end

% function [nl,nk] = GetNumberOfLine(fname)
%     % ... Öffne Datei mit Lesezugriff
%     fid = fopen(fname,'r');
%     % ... Einfach Zeilen Zählen bis ende erreicht
%     nl = 0;
%     stat = 1;
%     while stat
%         tline = fgetl(fid);
%         if tline == - 1
%             stat = 0;
%         else
%             if nl == 1
%                 k = strfind(tline,';');
%                 nk = length(k);
%             end
%             nl = nl + 1;
%         end
%     end
%     % ... Schließe Datei
%     fclose(fid);
% end
