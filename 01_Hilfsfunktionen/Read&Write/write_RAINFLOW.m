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
if size(P,1) == 6
    header = [' Mode  ',...
              ' Schwingspiel ',...
              '         P         ',...
              '         DLZ           ',...
              '         D             ',...
              '     Dsum       \n' ];
    formspec = '%3i;%10i;%20.7f;%20.7f;%20.7f;%20.7f\n';
    rowDLZ =4;
    rowP = 3;
    rowDsum = 6;
else 
    header = [' Schwingspiel ',...
              '            P          ',...
              '         DLZ           ',...
              '         D             ',...
              '     Dsum       \n' ];
    formspec = '%10i;%20.7f;%20.7f;%20.7f;%20.7f\n'; 
    rowDLZ =3;
    rowP = 2;
    rowDsum = 5;
end

% ... Schreibe Header
fprintf(fid,header);

% ... schreibe Daten
for i = 1:size(P,2)
    fprintf(fid,formspec,P(:,i));
end
fprintf(fid,'#HCM\n');

% -------------------------------------------------------------------------
% Ausgabe Kollektiv

% ... Herrausfiltern letzter Durchlauf
if ~isempty(P)
    ndl = ceil(P(rowDLZ,end));       % Anzahl Durchlaufe
    idxLast = P(rowDLZ,:) > ndl -1;  % Indices aller Beiträge des letzten DL
    % Schadenssumme vor letztem Durchlauf
    Dsum_vor_letztem_DL = max(P(rowDsum,~idxLast));
    % Schadensumme des letzten Durchlaufs
    Dsum_last = P(rowDsum,end) - Dsum_vor_letztem_DL;
    Pmax = max(P(rowP,idxLast));    % Maximalwert P im letzten DL
    
    % ... Klassieren letzter Durchlauf
    [H, H2] = klassieren(P(rowP,idxLast)./Pmax,32,0);
    
    
    % ... Schreibe Schadenssummen und Kolletivhöchstwert und Kollektiv
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'Anzahl simulierter Durchläufe: %i\n',ndl);
    fprintf(fid,'Schadensumme Durchäufe 1 - %i: %.15f\n',ndl-1,Dsum_vor_letztem_DL);
    fprintf(fid,'Schadensumme Durchlauf %i    : %.15f\n',ndl,Dsum_last);
    fprintf(fid,'\n');
    fprintf(fid,'Schädiungsparameterkollektiv in Durchlauf %i\n',ndl);
    fprintf(fid,'Pmax : %.5f\n',Pmax);
    fprintf(fid,'  Klasse    Häufigkeit    Summenhäufigkeit\n');
    for i = 1 : length(H)
        fprintf(fid,'%6i;%13i;%13i\n',length(H)-i+1,H2(i),H(i));
    end
    
    % ... Schließe Datei
    fclose(fid);
end
end % Ende Funktion