function write_RAINFLOW(jobname,namenszusatz,outpath,P,phi,psi)
% Schreibe Ergebnisse der Rainflowzählung und Schädigungsbewertung
% (ausrechnen von Schädigungsparametern) in eine Datei
%
% INPUT:
% jobname      - (str) Name der Rechnung
% namenszusatz - (str) zusätzlicher Name
% outpath   - (str) Pfad zu OutputFolder
% dmgmodel  - (int) Unterscheide Schädigungsmodelle
% P         - (struct) Ergebnisse der rainflowzählung (teilweise
%             abhängig von dmgmodel) Feldname = name des
%             Schädigungsparameters
% phi,psi   - (double) Aktuell betrachtete Ebene (in rad) 
%
% OUTPUT:
%  erzeugte Datei jobname_namenszusatz_PHI_PSI.hcm
%__________________________________________________________________________

% ... Umrechnen phi,psi in grad
phi = phi * 180/pi;
psi = psi * 180/pi;

% ... Dateiname
dateiname = [jobname,'_',namenszusatz,'_',num2str(phi),'_',num2str(psi),'.hcm'];
name = [outpath,'\',dateiname];

% ... Öffne Datei
fid = fopen(name,'w');

% ... def header und Format fürs Speichern von P
if size(P,1) == 4
    header = [' Mode  ',...
              ' Schwingspiel ',...
              '         P         ',...
              '         DLZ        \n' ];
    formspec = '%3i;%10i;%20.7f;%20.7f\n';      
else 
    header = [' Schwingspiel ',...
              '            P          ',...
              '         DLZ        \n' ];
    formspec = '%10i;%20.7f;%20.7f\n'; 

end

% ... Schreibe Header
fprintf(fid,header);

% ... schreibe Daten
for i = 1:size(P,2)
    fprintf(fid,formspec,P(:,i));
end

% ... Schließe Datei
fclose(fid);
end % Ende Funktion