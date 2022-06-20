function write_report(jobname,outpath...
                      ) 
% Schreibe Report der Rechnung in Datei. Enthält Übersicht über die
% getroffenten Einstellungen & Die Rechenergebnisse
%
% INPUT:
% jobname   - (str) Name der Outputdatei = jobname_SIGEPS.dat
% outpath   - (str) Pfad zu OutputFolder
%
% OUTPUT:
%  erzeugte Datei jobname.rpt
%
%__________________________________________________________________________

% ... Dateiname
dateiname = [jobname,'.cpl'];
name = [outpath,'\',dateiname];

% ... Öffne Datei
fid = fopen(name,'w');

% ... Schließe Datei
fclose(fid);

end % Ende Funktion



