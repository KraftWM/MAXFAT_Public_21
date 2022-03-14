function [SIG,EPS,DLZ] = read_sigepsfile(sigepsfile,ntens,ndi)
% Auslesen von binary Dateien
% gedacht fuer Ergebnisse der Materialmodelle
%
% INPUT:
% sigepsfile    - Dateiname
% ntens         - Anzahl Tensorkomponenten
% ndi           - Anzahl Nebendiagonalelemente
%
% OUTPUT:
% SIG           - Spannungen
% EPS           - Dehnungen
% DLZ           - Durchlaufzaehler
% _________________________________________________________________________


% oeffne Datei
fid = fopen(sigepsfile,'r'); 
% Unterscheide Spannungszustaende 
if ntens == 6 && ndi == 3 % 3D
    nin = 13; % DATA = fread(fid,[13,Inf],'double');
    nsig = 6;
    neps = 6;
elseif ntens == 3 && ndi == 2 % ESZ
    nin = 8; % DATA = fread(fid,[8,Inf],'double');
    nsig = 3;
    neps = 4;
elseif ntens == 2 && ndi == 1 % Sigma - Tau
    nin = 7; % DATA = fread(fid,[7,Inf],'double');
    nsig = 2;
    neps = 4;
elseif ntens == 1 && ndi == 1 % reiner Zug
    nin = 5; % DATA = fread(fid,[5,Inf],'double');
    nsig = 1;
    neps = 3;
else % Spannungszustand nicht erkannt
    msg = 'Spannungszustand nicht erkannt';
    error(msg)
end
% Lese Datei
DATA = fread(fid,[nin,Inf],'double');
% Schließe Datei
fclose(fid);
% Aufteilen
DLZ = DATA(1,:);
SIG = DATA(2:nsig+1,:);
EPS = DATA(nsig+2:1+nsig+neps,:);
end