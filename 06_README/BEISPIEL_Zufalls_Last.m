% ########################################################################
% # Beispiel: Berechnen einer zufällig generierten Lastfolge             #
% #                                                                      #
% # Werkstoff: Baustahl S355                                             #
% # Bauteil  : Abgesetzte Welle, Kerbradius 0.5mm                        #
% # Belastung: 90° Phasenverschiebung                                    #
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
jobname = 'Random_Load';                   % Name der Rechnung

if ~exist(outpath,'dir')
    mkdir(outpath)
end

%% Definition Materialmodell

E = 204000;                            % E-Modul 
nu = 0.3;                              % Querdehnzahl 
Rm = 541;                              % Zugfestigkeit 
% zyklisch stabil (Ramberg Osgood)
Ks = 839;% 1141;%                      % zyklische Steifigkeit 
ns = 0.138;%    0.193;%                % Verfestigungsexponent 
% Dehnungswoehlerlinie
sf = 762; ef = 0.415; b = -0.08; c = -0.556;
ND=5e5; 

%% Definitionen Bauteil

% Übertragungsfaktoren Kerbradius 0.5mm 
%     xx     yy     xy
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
eindkerb = 'Neuber Stern';
% Nichtproportionale Verfestigung
nppara = 'Riess';
% Materialmodell
material = 'KarimOhno';
% Anzahl der zu simulierenden Durchläufe durch die Lastfolge
ndl = 3;

%%  Definition Lastfolgen
% Fiktive Zeit
Sekunden=10;                                     % Dauer Signal
Fs=1200;                                         % Aufnahmefrequenz
ndp = fix(Sekunden*Fs);                          % Anzahl Datenpunkte
t=linspace(0,Sekunden,ndp);                      % fiktive Zeit
% zufällige Lastfolge der Nettospannung Zug 
Lamp = 100*rand(3,1);                            % Amplituden
phi = 2*pi*rand(3,1);                            % Phasenverschiebungen
Lmit = 50*rand(1)-25;                            % Mittelwert
omega = 2*pi*randi(Sekunden,3,1);                % (Kreis)Frequenzen
LZ = Lmit + sum(Lamp.*cos(omega.*t+phi),1);      % Lastfolge Zug
% zufällige Lastfolge der Nettospannung Torsion
Lamp = 100*rand(3,1);                            % Amplituden
phi = 2*pi*rand(3,1);                            % Phasenverschiebungen
Lmit = 50*rand(1)-25;                            % Mittelwert
omega = 2*pi*randi(Sekunden,3,1);                % (Kreis)Frequenzen
LT = Lmit + sum(Lamp.*cos(omega.*t+phi),1);      % Lastfolge Torsion
% Zusammenfassung
last = [LZ;LT];

%% Abschaetzen der Kennwerte nach FKM

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

% Absicherung Lastfolge
gamma_L1 = AbsicherungLastfolgeFKM(1,s_L,max(last(1,:)),P_L,P_A);
gamma_L2 = AbsicherungLastfolgeFKM(1,s_L,max(last(2,:)),P_L,P_A);
last = [gamma_L1; gamma_L2] .* last;
% Definition einer Instanz der Kerbsimulation
K = Kerbsimulation(jobname,outpath,...
    last,cf,ndl,Rm,E,nu,Ks,ns,...
    'material',material,...
    'eindkerb',eindkerb,...
    'nppara',nppara,...
    'Kp',Kp);
% Definition einer Instanz des Parameters PRAM
% 1. Pram aus Dehnungs-WL
% 2. Pram_stuetz anhand abgeschaetzter Stuetzstellen
Pram = PRAM(E,sf,ef,b,c,ND,Msig,'nst',n);
Pram_s = PRAM_stuetz(E,Pram_WS_stuetz,d1_ram,d2_ram,Msig,fram);
% Definition einer Instanz des Pramaters PFS
% 1. Pfs aus Dehnungs-WL
% 2. Pfs_stuetz anhand abgeschaetzter Stuetzstellen
Rp02s = Ks*0.002^ns;
sigF = 0.5 *(Rm+Rp02s);  % Fließgrenze als Mittelwert aus Rp02 und Rm
Pfs = PFS(E,nu,sigF,sf,ef,b,c,ND,'nst',n);
Pfs_var = PFS(E,nu,sigF,sf,ef,b,c,ND,'nst',n,'kfsopt','var','tf',sf/sqrt(3),'gf',ef*sqrt(3),'bg',b,'cg',c);
Pfs_s = PFS_stuetz(E,nu,Pfs.kfs,Pfs_WS_stuetz,d1_fs,d2_fs,fram);
% Definition einer Instanz des Pramaters PZ
% 1. Pz aus Dehnungs-WL
% 2. Pz_stuetz anhand abgeschaetzter Stuetzstellen
Pz = PZ2(E,nu,Rm,Ks,ns,sf,ef,b,c,ND,'nst',nst,'Gsig',Gsig,'Gtau',Gtau);
Pz_s = PZ2_stuetz(E,nu,sigF,Pz_WS_stuetz,d_z,fraj,...
    'Pz_WSD_stuetz',Pz_WSD_stuetz,'Gsig',Gsig,'Gtau',Gtau);
% Zusammenfassen aller Instanzen der definierten Schaedigungsparameter
% als cell array
DMGs = {Pram, Pram_s, Pfs, Pfs_var, Pfs_s, Pz, Pz_s};
%     DMGs = {Pram_s, Pfs_s, Pz_s};
% winkel = [phi_max, phi_min, dphi, psi_max, psi_min, dpsi]
winkel = [90 0 9 45 0 45];
% Funktion zum berechnen der Anrisslebensdauer
[DL,DLc,phic,psic] = schadensrechnung_kerb(...
    jobname,outpath,...% Name der Rechnung
    K,...              % Instanz der Kerbsimulation
    DMGs,...           % Cell array Instanzen Schädigungsparameter
    winkel,...         % Winkel kritische Ebenen Rechnung
    1,...              % Display ausgabe
    1,...              % File Ausgabe Ergebnisse kritische Ebene
    1,...              % Ausgabe Ergebnisse Rainflow Zählung
    0);                % Rainflow Zählung in allen Ebenen


% Ende der Rechnung

%% Plot der Ergebnisse

% Lastfolge
% Plot
figure, grid on, hold on
plot(LZ,'LineWidth',1.3);
plot(LT,'LineWidth',1.3);
ylabel('Last'),xlabel('fiktive Zeit'),title('zufällige Lastfolge')
legend('$L_{Zug}$','$L_{Torsion}$','Interpreter','latex','FontSize',14)