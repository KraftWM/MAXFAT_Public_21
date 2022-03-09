function [DLZ,SIG,EPS,EPSP] = read_SIGEPS(fname)
% Funktion ließt Spannungs und Dehnungswerte aus .sig datei ein
%
% INPUT:
% fname     - Dateiname
%
% OUTPUT:
% SIG       - Verlauf Spannungen
% EPS       - Verlauf Dehnungen
% EPSP      - Verlauf plastische Dehnungen
% DLZ       - Durchlaufzähler
%__________________________________________________________________________

% ... Öffne Datei mit Lesezugriff
fid = fopen(fname,'r');

% ... Anzahl Zeilen&Spalten in Datei
tline = fgetl(fid);
A = sscanf(tline,'%i;%i');
nrows = A(1);
ncols = A(2);

% Format
formspec = [];
for i = 1:ncols-1
    formspec = [formspec,'%f;'];
end
formspec = [formspec,'%f'];

% ... Interpretiere Header
tline = fgetl(fid);
ntens = length(strfind(tline,'S'));
plastopt = contains(tline,'P');

% ... Speicher für OUTPUT
DATA = fscanf(fid,formspec,[ncols,nrows]);
DLZ = DATA(1,:);
SIG = DATA(2:ntens+1,:);
if plastopt
    EPS = DATA(ntens+2:2*ntens+1,:)./100;
    EPSP = DATA(2*ntens+2:3*ntens+1,:)./100;
else
    EPS = DATA(ntens+2:2*ntens+2,:)./100;
    EPSP = NaN;
end
% ... Schließe Datei
fclose(fid);

end % Ende Funktion

