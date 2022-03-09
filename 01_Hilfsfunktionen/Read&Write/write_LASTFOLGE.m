function write_LASTFOLGE(jobname,outpath,lastfolge)
% Schreibe Spannungen und Dehnungen im Kerbkoordinatensystem in Datei
%
% INPUT:
% jobname   - (str) Name der Outputdatei = jobname_SIGEPS.dat
% outpath   - (str) Pfad zu OutputFolder
% Lastfolge - (double array) Lastfolge
%             lastfolge e R^(nkana x ndata)
%             
% 
% OUTPUT:
%  erzeugte Datei jobname_LASTFOLGE.lfo
%
%__________________________________________________________________________

% ... Dateiname
dateiname = [jobname,'.lfo'];
name = [outpath,'\',dateiname];

% ... Öffne Datei
fid = fopen(name,'w');

% ... anzahl daten
nkana = size(lastfolge,1);
ndata = size(lastfolge,2);

% ... Schreibe Header und Format
fprintf(fid,'Num Data: %i\nNum Kana: %i\n',ndata,nkana);
header = [];
formspec = [];
for i = 1 : nkana
    header = [header, '    Kanal ',num2str(i),'    '];
    formspec = [formspec,'%15.5f;'];
end
header = [header,'\n'];
formspec = [formspec,'\n'];

fprintf(fid,header);


% ... Schreibe Daten in Datei
for i = 1:ndata
    fprintf(fid,formspec,lastfolge(:,i));
end

% ... Schließe Datei
fclose(fid);
end % Ende Funktion