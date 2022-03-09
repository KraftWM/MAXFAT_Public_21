function [PHI,PSI,DL,zmin] = read_CRITPLANE(fname)
% Auslesen von .cpl Dateien
%
% INPUT:
% fname         - Dateiname
%
% OUTPUT:
% PHI           - Winkel
% PSI           - Winkel
% DL            - Durchläufe
% zmin          - Zeile kritische Ebene
% _________________________________________________________________________

% ... Anzahl (ununterbrochene) Zeilen in Datei
[nl] = GetNumberOfLine(fname);
% ... Öffne Datei mit Lesezugriff
fid = fopen(fname,'r');

% ... Anzahl Schädigungsparameter
tline = fgetl(fid);
np = length(strfind(tline,'('));

% ... Speicher für OUTPUT
PHI = zeros(nl-1,1);
PSI = zeros(nl-1,1);
DL  = zeros(nl-1,np);
% ... geringste Lebensdauern
DLmin = repelem(1e40,np);
zmin = repelem(0,np);



% ... Lese Datei
status = 1;
z = 1;
while status
    % ... Lese nächste Zeile
    tline = fgetl(fid);
    if tline == -1
        % ... Lesen abbrechen
        status = 0;
    elseif strcmp(tline,'')
        status = 0;
    else
        % ... lese Werte
        [PHI(z),PSI(z),DL(z,:),z] = readVal(tline,z,np);
        % Kleinste Werte
        for k = 1 : np
            if DL(z-1,k) < DLmin(k)
                DLmin(k) = DL(z-1,k);
                zmin(k) = z-1;
            end
        end
    end
end

% ... Schließe Datei
fclose(fid);

end % Ende Funktion


% ----------------------------------------------------------------------- %
% Hilfsfunktion
% ----------------------------------------------------------------------- %
function [phi,psi,dl,z] = readVal(tline,z,np)
    k = strfind(tline,';');
    phi = str2double(tline(1:k(1)-1));
    psi = str2double(tline(k(1)+1:k(2)-1));
    dl  = zeros(1,np);
    r = 1;
    for i = 2 :  length(k) - 1        
        dl(r)  = str2double(tline(k(i)+1:k(i+1)-1));
        r = r + 1;
    end
    dl(r)  = str2double(tline(k(i+1)+1:end));
    z = z + 1;
end

function [nl] = GetNumberOfLine(fname)
    % ... Öffne Datei mit Lesezugriff
    fid = fopen(fname,'r');
    % ... Einfach Zeilen Zählen bis ende erreicht
    nl = 0;
    stat = 1;
    while stat
        tline = fgetl(fid);
        if tline == - 1
            stat = 0;
        elseif strcmp(tline,'')
            stat = 0;
        else
            nl = nl + 1;
        end
    end
    % ... Schließe Datei
    fclose(fid);
end