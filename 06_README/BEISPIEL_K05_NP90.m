% ########################################################################
% # Beispiel: Berechnen einer Bauteilfließkurve                          #
% #                                                                      #
% # Werkstoff: Baustahl S355                                             #
% # Bauteil  : Abgesetzte Welle, Kerbradius 0.5mm                        #
% # Belastung: 90° Phasenverschiebung, R = -1                            #
% #            Lastkanal 1 - Zug/Druck                                   #
% #            Lastkanal 2 - Torsion                                     #
% ########################################################################

% Clear Matlab
clear 
clc
% Füge alle Funktionen Matlab Pfad hinzu
% run ../startup

%% Definitionen Rechnung

outpath = '../00_Temp';                % Name des Output Folders
jobpart = 'K05_NP_';                   % Name der Rechnung

if ~exist(outpath,'dir')
    mkdir(outpath)
end

%% Definition Materialmodell

E = 204000;                            % E-Modul 
nu = 0.3;                              % Querdehnzahl 
Rm = 541;                              % Zugfestigkeit 
% zyklisch stabil (Ramberg Osgood)
Ks = 839;% 1141;%                      % zyklische Steifigkeit 
ns = 0.138;%    0.193;%                % verfestigungsexponent 
% Dehnungswoehlerlinie
sf = 762; ef = 0.415; b = -0.08; c = -0.556;
ND=5e5; 

%% Definitionen Bauteil

% Übertragungsfaktoren Kerbradius 0.5mm 
%     xx     yy     yy
cf = [3.351 0.8351 0.0;...            % Nettospannung Zug
      0.0   0.0    1.963]';           % Nettospannung Torsion
Kp = 2.7;                             % Traglastformzahl
Asig = 25;                            % Hochbeanspruchte Oberflaeche
Gsig = 2/0.5;                         % Gradient Normalspannung
Gtau = 1/0.5 + 2/25;                  % Gradient Schubspannung
Gmises = sqrt(Gsig^2+3*Gtau^2);              
Rz = 0;                               % Rauheit (unberücksichtigt)

%% Absicherung 

% Ausfallwahrscheinlichkeit (Tabelle 2.4 Richtlinie Nichtlinear)
P_A = 0.5;
% Auftretenswahrscheinlichkeit der Lastfolge
P_L = 0.5;
% Standartabweichung Absicherung der Lastfolge
s_L = 0;


%% Definition Kerbsimulation

% Bestimmen der Bauteilfließkurve
eindkerb = 'Seeger Beste';
% Nichtproportionale Verfestigung
nppara = 'Gaier';
% Anzahl der zu simulierenden Durchläufe durch die Lastfolge
ndl = 10;

%%  Definition Lastfolgen

% Anzahl Datenpunkte pro Zyklus
ndp = 64;
% Lastfolge der Nettospannung Zug mit Amplitude 1 
LZ = lastgenerator(1,1,90,1,0,ndp);  % = cos(2*pi*t);              
% Lastfolge der Nettospannung Torsion mit Amplitude 1
LT = lastgenerator(1,1,0,1,0,ndp);  % = sin(2*pi*t);
% Zusammenfassung
L = [LZ;LT];
% Amplituden der Nettospannungen aus den Versuchen
%         Zug/Druck   Torsion
AMP = [ ...
        203.7183272	181.0829575;...
        142.602829	126.7580702;...
        81.48733086	72.43318299;...
        101.8591636	90.54147874;...
        122.2309963	108.6497745;...
        162.9746617	144.866366;...
        224.0901599	199.1912532;...
        183.3464944	162.9746617...
      ]';

%% Abschaetzen Kennwerte nach FKM

% Werkstoffgruppe
wsgruppe = 'Stahl';
% Rauheitsfaktor
KPR = RauheitseinflussFKM(wsgruppe,Rz,Rm);
% Stuetzwirkung Richlinie Nichtlinear
[n,nst,nbm] = StuetzZifferFKM(Asig,Gmises,Rm,wsgruppe);
% Pram - Woehlerlinie
Msig = MittelspannungseinflussFKM(wsgruppe,Rm);
[Pram_WS_stuetz,Pram_WSD_stuetz,d1_ram,d2_ram,f0025] = ...
                      Pram_Woehlerlinie_FKM(Rm,wsgruppe);
% Absicherung Bauteilwoehlerlinie
fram = AbsichernBauteilWoehlerlinie_Pram(P_A,n,KPR,f0025);
% Pfs - Woehlerlinie (vorlaufig)
[Pfs_WS_stuetz,Pfs_WSD_stuetz,d1_fs,d2_fs] = ...
                      Pfs_Woehlerlinie(Rm,wsgruppe);
% Pz - Woehlerlinie (vorlaufig aus Praj Woehlerlinie)
[Pz_WS_stuetz,Pz_WSD_stuetz,d_z,f0025] = ...
                      Pz_Woehlerlinie(Rm,wsgruppe);
% Absicherung Bauteilwoehlerlinie
fraj = AbsichernBauteilWoehlerlinie_Praj(P_A,n,KPR,f0025);

%% Rechnung 

% Zum berechnen der Bauteilwöhlerlinie werden alle Amplituden aus der
% Versuchsreihe nachgerechnet. Dazu wird in jedem Schleifendurchlauf die
% Lastfolge skaliert. Da alle Rechnungen unabhängig von einander sind,
% können diese auch parallel ausgeführt werden.

% parfor i = 1:size(AMP,2) % Parallele Rechnung, ! optdisplay auf 0 setzten

for i = 1:size(AMP,2) % Sequenzielle Rechnung    
    % Amplitude aktueller Versuch
    amp = AMP(:,i);
    last = amp .* L;
    % Absicherung Lastfolge
    gamma_L1 = AbsicherungLastfolgeFKM(1,s_L,max(last(1,:)),P_L,P_A);
    gamma_L2 = AbsicherungLastfolgeFKM(1,s_L,max(last(2,:)),P_L,P_A);
    last = [gamma_L1; gamma_L2] .* last;
    % Jobname 
    jobname = [jobpart,num2str(amp(1))];
    % Definition einer Instanz der Kerbsimulation
    K = Kerbsimulation(jobname,outpath,...
                last,cf,ndl,Rm,E,nu,Ks,ns,...
                'eindkerb',eindkerb,...
                'nppara',nppara,...
                'Kp',Kp);
    % Definition einer Instanz des Parameters PRAM
    % 1. Pram aus DehnungsWL
    % 2. Pram anhand abgeschaetzter Stuetzstellen
    Pram = PRAM(E,sf,ef,b,c,ND,Msig,'nst',n);
    Pram_s = PRAM_stuetz(E,Pram_WS_stuetz,d1_ram,d2_ram,Msig,fram);
    % Definition einer Instanz des Pramaters PFS
    % 1. Pfs aus DehnungsWL
    % 2. Pfs anhand abgeschaetzter Stuetzstellen
    Rp02s = Ks*0.002^ns;
    sigF = 0.5 *(Rm+Rp02s);  % Fließgrenze als Mittelwert aus Rp02 und Rm
    Pfs = PFS(E,nu,sigF,sf,ef,b,c,ND,'nst',n);
    Pfs_s = PFS_stuetz(E,nu,Pfs.kfs,Pfs_WS_stuetz,d1_fs,d2_fs,fram);
    % Definition einer Instanz des Pramaters PZ
    % 1. Pz aus DehnungsWL
    % 2. Pz anhand abgeschaetzter Stuetzstellen
    Pz = PZ2(E,nu,Rm,Ks,ns,sf,ef,b,c,ND,'nst',nst,'Gsig',Gsig,'Gtau',Gtau);
    Pz_s = PZ2_stuetz(E,nu,sigF,Pz_WS_stuetz,d_z,fraj,...
           'Pz_WSD_stuetz',Pz_WSD_stuetz,'Gsig',Gsig,'Gtau',Gtau);
    % Zusammenfassen aller Instanzen der definierten Schaedigungsparameter
    % als cell array
    DMGs = {Pram, Pram_s, Pfs, Pfs_s, Pz, Pz_s};
%     DMGs = {Pram_s, Pfs_s, Pz_s};
    % winkel = [phi_max, phi_min, dphi, psi_max, psi_min, dpsi]
    winkel = [90 0 9 45 0 45];
    % Funktion zum berechnen der Anrisslebensdauer
    [DL,DLc,phic,psic] = schadensrechnung_kerb(...
              jobname,outpath,...   % Name der Rechnung
              K,...           % Instanz der Kerbsimulation
              DMGs,...        % Cell array Instanzen Schädigungsparameter
              winkel,...      % Winkel kritische Ebenen Rechnung
              1,...           % Display ausgabe
              1,...           % File Ausgabe Ergebnisse kritische Ebene
              1,...           % Ausgabe Ergebnisse Rainflow Zählung
              0);             % Rainflow Zählung in allen Ebenen
end

% Ende der Rechnung

%% Auswerten der Rechnung 
% Lese alle Ergebnisse
list = dir([outpath,'\*',jobpart,'*.cpl']);
Results = NaN(length(list),19);
for i = 1 :length(list)
    % Lese Ergebnisse
    [PHI,PSI,DL,z] = read_CRITPLANE([list(i).folder,'/',list(i).name]);
    % Amplitude
    ampstr = list(i).name;
    ampstr = erase(ampstr,jobpart);
    ampstr = erase(ampstr,'.cpl');
    amp = str2double(ampstr);
    % Speichern Ergebnisse
    Results(i,1) = amp;
    for j = 1 : length(z)
        Results(i,(j-1)*3+2:j*3+1) = [PHI(z(j)),PSI(z(j)),DL(z(j),j)];
    end
end
[~,I] = sort(Results(:,1),1,"ascend");
Results = Results(I,:);

% Versuchsergebnisse
Mess = [ ...
203.7183272	1600;...
142.602829	5800;...
81.48733086	NaN;...
101.8591636	32500;...
122.2309963	21500;...
162.9746617	4800;...
224.0901599	1680;...
183.3464944	3600];

% Plote Ergebnisse
mycol = lines;
% PRAM
myplot(4,7,'RAM',Results,Mess,mycol(1,:));
% PFS
myplot(10,13,'FS',Results,Mess,mycol(2,:));
% PZ
myplot(16,19,'Z',Results,Mess,'k');


% Dummy Funktion zum ploten
function ax = myplot(col,col_s,namePara,Results,Mess,linecol)
lw = 2;
colorAnriss = [0.00, 0.549, 0.3098];
material = 'S355';
rho = '$\rho = 0.5$ mm';
last = 'Nichtprop. $R=-1$';
ax = plot_BWL(0,Mess(:,2),Mess(:,1),1,'','Messung','','Marker','o','Color',colorAnriss,'LineWidth',lw,'LineStyle','none');
plot_BWL(ax,Results(:,col),Results(:,1),0,'',['$P_{',namePara,'}$'],'','Color',linecol,'LineWidth',lw,'LineStyle','-');
plot_BWL(ax,Results(:,col_s),Results(:,1),0,'',['$P_{',namePara,',stuetz}$'],'','Color',linecol,'LineWidth',lw,'LineStyle','--');
title(['$P_{',namePara,'}$ Bauteilwoehlerlinie'])
% Textbox
dimHead = [0.145,0.25,0.3,0.06];
dimBody = [0.145,0.165,0.3,0.085];
ampstr = {'\textbf{Belastung}:',last};
annotation('textbox',dimHead,'String',['\textbf{',material,'} ',rho],'FitBoxToText','off',"Interpreter","latex","FontName",'Arial','FontSize',13);
annotation('textbox',dimBody,'String',ampstr,'FitBoxToText','off',"Interpreter","latex","FontName",'Arial');
end