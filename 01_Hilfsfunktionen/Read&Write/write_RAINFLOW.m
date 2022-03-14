function write_RAINFLOW(jobname,namenszusatz,outpath,P,phi,psi)
% Schreibe Ergebnisse der Rainflowzaehlung und Schaedigungsbewertung
% (ausrechnen von Schaedigungsparametern) in eine *.hcm Datei
%
% INPUT:
% jobname      - (str) Name der Rechnung
% namenszusatz - (str) zusaetzlicher Name
% outpath   - (str) Pfad zu OutputFolder
% dmgmodel  - (int) Unterscheide Schaedigungsmodelle
% P         - (struct) Ergebnisse der rainflowzaehlung (teilweise
%             abhaengig von dmgmodel) Feldname = name des
%             Schaedigungsparameters
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

% ... oeffne Datei
fid = fopen(name,'w');

% ... def header und Format fuers Speichern von P
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