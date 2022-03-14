function [RESULTS,DLc,phic,psic,Totaltime] = schadensrechnung_mat(jobname,outpath,...
                                    material,para,M,load,inkflag,ntens,ndi,...
                                    ndl,...
                                    DMGs,...
                                    winkel,...
                                    optdisplay,...
                                    optcritplane,...
                                    optrainflow,...
                                    optallhcm)
% Hauptfunktion zum ausfuehren der Schaedigungsrechnung fuer Material
%
% BERECHNUNGSABLAUF:
%
% Intergration Materialmodell    - Berechne lokale Beanspruchung
%
% Kritische Ebenen Schleife      - Transformation lokaler Spannungen und
%                                  Dehnungen in verschiede Schnittebenen 
%                                  und Schaedigungsrechnng in den
%                                  Schnittebenen
%
% Rainflow (HCM)                 - Zyklenzaehlen in den verschiedenen Ebenen
%
% Schaedigungsrechnung            - Sobald eine Hystere schließt wird direkt
%                                  ein Schaedigungsparameter berechnet
%
% -------------------------------------------------------------------------
%
% INPUT:
% 1. Dateiverwaltung
% jobname          - (str) name der rechnung
% oupath           - (str) pfad indem output dateien gespeichert werden
%
% 3. Integration Materialmodell
% material         - (str) name des Materialmodells 
%                      'Chaboche','KarimOhno','OhnoWang'
% para             - (double array) Parameter des Materialmodells
%  M               - (int) Anzahl Backstresstensoren
% load             - (double array) Last-Zeitreihe
% inkflag          - (int) 0 load enthaelt Spannungen, 1 load enthaelt Dehnungen
% ntens,ndi        - (int) Anzahl Tensorkomponenten zum Unterscheiden der
%                     Spannnungszustaende
% ndl              - (int) Anzahl zu simulierender Durchlaeufe
%
% 4. Schaedigung
% DMGs             - (cell array mit Objekten von Schaedigungsparametern)
%
% 5. kritische Ebene 
% winkel           - (double array), Winkel fuer kritische Ebenen Schleife
%                    in Grad
%                    phimax = winkel(1);
%                    phimin = winkel(2);
%                    dphi   = winkel(3);
%                    psimax = winkel(4);
%                    psimin = winkel(5);
%                    dpsi   = winkel(6);
%
% 6. Outputoptionen
% optdisplay       - Anzeige auf dem Display
% optcritplane     - Ausgabe der Ergebnisse ALLER Ebenen in
%                               Datei
% optrainflow      - Schreibe Ergebnisse der Rainflowzaehlung in
%                               eine Datei
% optallhcm        - Schreibe Ergebnisse der Rainflowzaehlung in
%                               ALLEN Ebenen in eine Datei
%
% OUTPUT:
% RESULTS         - (double array) Ergebnisse der Rechnung
%                   1. Spalte phi der Ebene 
%                   2. Spalte psi der Ebene
%                   3.-... Spalte Durchlaeufe bis Anriss fuer verschiedene 
%                          Schaedigungsmodelle
% DLc             - Durchlaufe in kritischer Ebene
% phic,psic       - Winkel kritische Ebene
% Totaltime       - Dauer der Gesamten Rechnung (nur wenn Displayausgabe an
%                   ist, sonst Totaltime = 0)
%__________________________________________________________________________


% Displayausgabe
if optdisplay
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('Jobname           :%s\n',jobname);
    fprintf('Materialmodell    :%s\n',material);
    fprintf('Backstresstensoren:%i\n',M);
    fprintf('Anzahl Druchlaeufe :%i\n',ndl);
    fprintf('Anzahl Inkemente  :%i\n',length(load));
    fprintf('Starte Integration Materialmodell, ');
    tic;
end
Totaltime = 0;       % Runtime

% Integration Materialmodell
[sigepsfile,ndata]  = solve_matmodel_file(jobname,outpath,...
                         material,ntens,ndi,para,M,inkflag,load,ndl);

% Displayausgabe
if optdisplay
    time = toc;
    Totaltime = Totaltime + time;
    fprintf('Ende, Dauer: %.5f\n',time);
end

%% Kritische Ebenen Rechnung
% Winkel fuer kritische Ebene
dphi = winkel(3); phimin = winkel(2); phimax = winkel(1);
dpsi = winkel(6); psimin = winkel(5); psimax = winkel(4);
if optdisplay
    numpsi = ceil((psimax-psimin)/dpsi)+1;
    numphi = ceil((phimax-phimin)/dphi)+1;
    numwinkel = numpsi * numphi;
    fprintf('--------------------------------------------------------------------------\n')
    fprintf('Starte kritische Ebenen Rechnung\n');
    fprintf('%5s%7.3f%9s%7.3f%9s%7.3f\n','dphi:',dphi,', phimin:',phimin,', phimax:',phimax);
    fprintf('%5s%7.3f%9s%7.3f%9s%7.3f\n','dpsi:',dpsi,', psimin:',psimin,', psimax:',psimax);
    fprintf('Anzahl zu berechnender Ebenen:%i\n',numwinkel);
    tic;
end

% kritische Ebenen Rechnung - Rainflow fuer Parameter Zusammen druchfuehren
[phic,psic,DLc,RESULTS] = criticalplaneV6(...
                     sigepsfile,ndata,ntens,ndi,...
                     dphi,phimax,phimin,dpsi,psimax,psimin,...
                     DMGs,...
                     optdisplay,...
                     optcritplane,...
                     optrainflow,...
                     optallhcm,...
                     jobname,outpath);
                 
if optdisplay
    time = toc;
    Totaltime = Totaltime + time;
    fprintf('Ende, Dauer krtische Ebene: %.5f\n',time);
    fprintf('      Dauer gesamt        : %.5f\n',Totaltime);
    fprintf('kritische Ebene:\n');
    fprintf('Parameter     phi     psi      DL\n');
    for i = 1 : length(DMGs)
        fprintf('%8s:%8.3f%8.3f%15.3f\n',DMGs{i}.Name,phic(i),psic(i),DLc(i));
    end
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
end


end % Ende Funktion