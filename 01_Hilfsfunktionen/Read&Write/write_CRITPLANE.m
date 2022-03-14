function write_CRITPLANE(jobname,outpath,dmgmodel,D) 
% Schreibe Ergebnisse der kritischen Ebenen Rechnung in Datei
%
% INPUT:
% jobname   - (str) Name der Outputdatei = jobname_SIGEPS.dat
% outpath   - (str) Pfad zu OutputFolder
% dmgmodel  - (cell array) objekte der Schaedigungsparameter
%  D        - (double array) 
%            D(:,1) = Winkel phi
%            D(:,2) = Winkel psi
%            D(:,3) = Durchlaeufe
%
% OUTPUT:
%  erzeugte Datei jobname.cpl
%
%__________________________________________________________________________

% ... Anzahl der Modelle
numdmgs = length(dmgmodel);

% ... Dateiname
dateiname = [jobname,'.cpl'];
name = [outpath,'\',dateiname];

% ... Oeffne Datei
fid = fopen(name,'w');

% ... Schreibe Header & Format
fprintf(fid,['        phi  ',...
             '        psi  ']);
formspec = '%12.3f;%12.3f';
for i = 1:numdmgs
    sizename = length(dmgmodel{i}.Name);
    fprintf(fid,['       DL(',dmgmodel{i}.Name,')']) ;
    for j = 1:19-sizename
        fprintf(fid,' ');
    end
    formspec = [formspec,';%30.7f'];
end
formspec = [formspec,'\n'];
fprintf(fid,'\n');

% ... schreibe Daten
for i = 1:size(D,1)
    fprintf(fid,formspec,D(i,:));
end

% ... Schreibe kritische Ebene
fprintf(fid,'\n\n\n');
fprintf(fid,'Kritische Ebene:\n');
fprintf(fid,'Parameter     phi     psi      DL\n');
for i = 1:numdmgs
    [DLcrit,I] = min(D(:,2+i));
    phic = D(I,1);
    psic = D(I,2); 
    fprintf(fid,'%8s:%8.3f%8.3f%15.3f\n',dmgmodel{i}.Name,phic,psic,DLcrit);
end 

% ... Schließe Datei
fclose(fid);

end % Ende Funktion



