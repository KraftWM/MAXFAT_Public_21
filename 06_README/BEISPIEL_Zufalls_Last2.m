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
% run startup

%% Definitionen Rechnung

outpath = '00_Temp';                       % Name des Output Folders
jobname = 'Random_Load';                   % Name der Rechnung

if ~exist(outpath,'dir')
    mkdir(outpath)
end

%% Definition Materialmodell

% % Aus Versuchen
E = 204000;                            % E-Modul 
nu = 0.3;                              % Querdehnzahl 
Rm = 541;                              % Zugfestigkeit 
% zyklisch stabil (Ramberg Osgood)
Ks = 839;                              % zyklische Steifigkeit 
ns = 0.138;                            % Verfestigungsexponent 
% Dehnungswoehlerlinie
sf = 762; ef = 0.415; b = -0.08; c = -0.556;
ND=5e5; 

% Abgeschätzt
% wsgruppe = 'Stahl';
% Rm = 541;
% [E,nu] = StatischFKM(wsgruppe);                                            % Elastisches Materialverhalten
% [Ks,ns] = RambergOsgoodFKM(wsgruppe,Rm);                                   % zyklisches Materialverhalten
% [sf,ef,b,c] = DehnungsWoehlerlinie_Waechter(wsgruppe,Rm);                  % DehnungsWL
% ND = 5e5;                                                                  % Dauerfestigkeit

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
ndl = 2;

%%  Definition Lastfolgen
load('GAUSS4_Filter03_NP.mat')

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
gamma_L1 = AbsicherungLastfolgeFKM(1,s_L,max(L(1,:)),P_L,P_A);
gamma_L2 = AbsicherungLastfolgeFKM(1,s_L,max(L(2,:)),P_L,P_A);
L = [200*gamma_L1; 200*gamma_L2] .* L;
% Definition einer Instanz der Kerbsimulation
K = Kerbsimulation(jobname,outpath,...
    L,cf,ndl,Rm,E,nu,Ks,ns,...
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
% 2. Pfs aus Dehnungs-WL kfs = kfs(N)
% 2. Pfs_stuetz anhand abgeschaetzter Stuetzstellen
Rp02s = Ks*0.002^ns;
sigF = 0.5 *(Rm+Rp02s);  % Fließgrenze als Mittelwert aus Rp02 und Rm
[tf,gf,bg,cg] = GleitungsWL(sf,ef,b,c,1); % Abschätzen GLeitungs-WL
Pfs = PFS(E,nu,sigF,sf,ef,b,c,ND,'nst',n);
Pfs_var = PFS(E,nu,sigF,sf,ef,b,c,ND,'nst',n,'kfsopt','var','tf',tf,'gf',gf,'bg',bg,'cg',cg);
Pfs_s = PFS_stuetz(E,nu,Pfs.kfs,Pfs_WS_stuetz,d1_fs,d2_fs,fram);
% Definition einer Instanz des Pramaters PZ
% 1. Pz aus Dehnungs-WL
% 2. Pz_stuetz anhand abgeschaetzter Stuetzstellen
Pz = PZ2(E,nu,Rm,Ks,ns,sf,ef,b,c,ND,'nst',nst,'Gsig',Gsig,'Gtau',Gtau);
Pz_s = PZ2_stuetz(E,nu,sigF,Pz_WS_stuetz,d_z,fraj,...
    'Pz_WSD_stuetz',Pz_WSD_stuetz,'Gsig',Gsig,'Gtau',Gtau);


% Zusammenfassen aller Instanzen der definierten Schaedigungsparameter
% als cell array
DMGs = {Pram, Pz};
% winkel = [phi_max, phi_min, dphi, psi_max, psi_min, dpsi]
winkel = [90 0 9 45 0 45];
% Funktion zum berechnen der Anrisslebensdauer
[DL,DLc,phic,psic] = schadensrechnung_kerb(...
    jobname,outpath,...% Name der Rechnung
    K,...              % Instanz der Kerbsimulation
    DMGs,...           % Cell array Instanzen Schädigungsparameter
    winkel,...         % Winkel kritische Ebenen Rechnung
    1,...              % Display Ausgabe
    1,...              % File Ausgabe Ergebnisse kritische Ebene
    1,...              % Ausgabe Ergebnisse Rainflow Zählung
    0);                % Rainflow Zählung in allen Ebenen


% Ende der Rechnung

%%
% [PHI,PSI,DL2,zmin] = read_CRITPLANE([outpath,'/',jobname,'.cpl']);
% plot_CRITPLANE(PHI,PSI,DL2(:,6),1,'$P_Z$')