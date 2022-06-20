function [ZVARneu] = OhnoWangImp(ntens, ndi, ink, ZVAR, ink_flag, para)
% Ohno Wang Plastizitätsmodell
%
% Aktuelle Implementierung:
% 3D - Dehnungssteuerung
% 
% TODO:
% 3D  - Spannungssteuerung
% ESZ - Dehnungssteuerung
% ESZ - Spannungssteuerung
% SubStepping Algo
% LineSearch Algo
% Ableitung des Residuums nicht numerisch sondern analystisch berechnen
%
% QUELLE MODELL:
% Aus Ohno et al. 1993 KINEMATIC HARDENING RULES WITH CRITICAL
% STATE OF DYNAMIC RECOVERY, PART I
%
% QUELLE Integration:
% Seifert,Schenk,Schmidt 2007 - Efficient and modular algorithms in 
% modeling finite inelastic deformations: Objective integration, parameter 
% identification and sub-stepping techniques
% 
% Seifert,Schmidt 2008 - Line-search methods in general return mapping
% algorithms with application to porous plasticity
% 
% Computational inelasticity
%
%   INPUT:
%         ntens     -> Anzahl Tensorkomponten (in Spannungen)
%         ndi       -> Anzahl Diagonalelemente, mehr kkommentare
%         ink       -> Belastungsinkrement
%         ZVAR      -> Zustandsvariablen [eps;epsp;alphai;p] bein spansteu
%                                        [sig;epsp;alphai;p] bein dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         para      -> Materialparameter des Modells
%                      [E,nu,c_i,r_i,chi_i,r0]
%
%
%    OUTPUT:
%        ZVARneu -> neue zustandsvariablen nach Lastinkrement
%        CEP     -> Tangentiale Steifigkeit
%    
%    SPANNUNGSZUSTÄNDE:
%      3D  - ntens=6, ndi=3
%     ESZ  - ntens=3, ndi=2
%    
%    NOTATIONEN (Beispiel ESZ):
%    Spannungen                 Verzerrungen
%
%           sig_11                     eps_11
%    sig =  sig_22               eps = eps_22
%           sig_12                     2eps_12
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2021 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: Januar 2021                                                  |
%  ----------------------------------------------------------------------

%--------------------------------------------------------------------------
%                   Unterscheide Spannungszustände
%--------------------------------------------------------------------------
% ... Elastische Parameter
E = para(1);
nu = para(2);

% ... 3D
if ntens == 6 && ndi == 3
    % ... Setze Abbildungen
    C = elast_steifigkeit(E,nu,ntens,ndi);
    D = elast_nachgiebigkeit(E,nu,ntens,ndi);
    [P,P2] = set_maps(ntens,ndi);
    
    % ... gehe in die Materialfunktion
    if ink_flag == 0 % Spannungssteuerung
    else % Dehnungssteuerung
        ZVARneu = OhnoWang3D_EPS(ZVAR,ink,para,C,P,P2);
    end
% ... ESZ    
elseif ntens == 3 && ndi == 2
    % ... Setze Abbildungen
    
    % ... gehe in die Materialfunktion
    
% ... Falsche Eingabe abfangen    
else
    error('Angegebener Spannungszustand in OhnoWang Modell ist falsch');
end

end % Ende Hauptfunktion



% ----------------------------------------------------------------------- %
% Integration 3D Spannungszustand                                         %
% ----------------------------------------------------------------------- %
function ZVAR1 = OhnoWang3D_EPS(ZVAR0,dEPS,para,...
                                C,P,P2)
% Integration 3D Spannungszustand
% INPUT:
%      ZVAR0  - Zustandsvariablen zum Start des Inkrements
%      dEPS   - Dehnungsinkrement
%      para   - Modellparameter
%      C      - Elastizitätstensor
%      P      - Abbildung Spannung auf Spannungsdeviator
%      P2     - Skalarprodukt Spannungsdeviatoren & Verdopplung Von
%               Nebendiagonalelementen bei Dehnungstensoren
% OUTPUT:
%      ZVAR1  - Zustandsvariablen am Ende des Inkrements
% 
% ZVAR = [sig;epsp;alphai;p]
% para = [E,nu,c_i,r_i,chi_i,r0]
% -------------------------------------------------------------------------

% Ein paar Definitionen ---------------------------------------------------
ntens = 6;                                                                 % Tensorkomponenten
M = (length(para)-3)/3;                                                    % Anzahl TBST
Y = para(end);                                                             % startradius fliessfläche
c = para(3:2+M);                                                           % Verfestigungsgeschwindigkeit
r = para(3+M:2+2*M);                                                       % Grenzradien Backstress
chi = para(3+2*M:2+3*M);                                                   % Ratchetting Exponent

% Zustandsvariablen -------------------------------------------------------
SIG0 = ZVAR0(1:ntens);                                                     % Spannungen
EPSP0 = ZVAR0(ntens+1:2*ntens);                                            % plastische Dehnungen
ALPHA0 = reshape(ZVAR0(2*ntens+1:(M+2)*ntens),ntens,M);                    % Backstresstensoren
p0 = ZVAR0((M+2)*ntens+1);                                                 % plastische Bogenlänge

% Elastischer Trial Strep -------------------------------------------------
SIG1 = SIG0+C*dEPS;
B1 = P * SIG1 - sum(ALPHA0,2);
F = MisesFF(B1,Y,P2);

% ... Trial Schritt wird angenommen
if F <= 0
    ZVAR1 = [SIG1; ZVAR0(ntens+1:end)];
    return;
end

% Trial Step wurde abgelehnt - Radial Return ------------------------------
% ... Konstanten für iteration
maxiter = 25;                                                              % maximale Iterationen
TOL = 1e-8;                                                                % Toleranz für Abbruchbedingung
iter = 0;                                                                  % Zähler durchgeführte Iterationen                                                                  
theta = 1;                                                                 % Faktor Liniensuche

% ... Initialisieren (Startwerte = Lösungen Trial Step)
EPSP1 = EPSP0;
ALPHA1 = ALPHA0;
dp1 = 0;
R = Residuum(dEPS,...                                                                         % Residuum
             SIG0,EPSP0,ALPHA0,...
             SIG1,EPSP1,ALPHA1,dp1,...
             ntens,C,P,P2,...
             Y,M,c,r,chi); 

% ... Schleife Newton Verfahren         
while abs(F) > TOL
    % ... Ableitung
    dRdY = VorwaertsDifferenzen(R,dEPS,...                                                                         % Residuum
                                SIG0,EPSP0,ALPHA0,...
                                SIG1,EPSP1,ALPHA1,dp1,...
                                ntens,C,P,P2,...
                                Y,M,c,r,chi,ntens*(M+2)+1);
    % ... Newtonschritt
    ZVAR1 = [SIG1; EPSP1; reshape(ALPHA1,M*ntens,1); dp1];
    ZVAR1 = ZVAR1 - theta * (dRdY\R);
    SIG1 = ZVAR1(1:ntens);                                                     
    EPSP1 = ZVAR1(ntens+1:2*ntens);                                            
    ALPHA1 = reshape(ZVAR1(2*ntens+1:(M+2)*ntens),ntens,M);    
    dp1 = ZVAR1((M+2)*ntens+1);
    % ... Prüfe Konvergenz
    R = Residuum(dEPS,...                                                                         % Residuum
             SIG0,EPSP0,ALPHA0,...
             SIG1,EPSP1,ALPHA1,dp1,...
             ntens,C,P,P2,...
             Y,M,c,r,chi); 
    F = R(ntens*(M+2)+1);
    % ... Inkrement Iterationszähler
    iter = iter + 1;
    % ... Abfangen Endlosschleife
    if iter > maxiter 
        disp('Scheiße');
        break;
    end
end % Ende Radial Return
% ... Anpassen plastische Bogenlänge
ZVAR1((M+2)*ntens+1) = ZVAR1((M+2)*ntens+1) + p0;

end % Ende Radial Return OhnoWang 3D





% -------------------------------------------------------------------------
% Auswerten der Fließfunktion nach von Mises
% -------------------------------------------------------------------------
function F2 = MisesFF(B,Y,P2)
% INPUT:
%  B      - effektive Spannung
%  Y      - Radius
%  P2     - Skalarprodukt deviatoren
% OUTPUT:
%  F2     - Fließfunktion zum Quadrat
%--------------------------------------------------------------------------                                   
F2 = 3/2 * B'*P2*B-Y^2;                                                           % Auswertung Fließfläche

end % Ende Mises Fließfunktion

% -------------------------------------------------------------------------
% Aufstellen des Flussvektors
% -------------------------------------------------------------------------
function N = FlussVektor(n,B,ALPHA,C,P2,Y,M,c,r,chi,sizeZ,ntens)
% INPUT:
%   n     - Normale an die Fließfläche
%   B     - effektive Spannung
%  ALPHA  - Backstresstensoren
%  C      - Elastizitätstensor
%  P2     - Skalarprodukt deviatoren
%  Y      - Radius Fließfläche
%  M      - Anzahl Tensorkomponenten
% c,r,chi - Parameter OW Modell
% sizeZ   - Anzahl Zustandsvariablen
% ntens   - Anzahl Tensorkomponenten
%
% OUTPUT:
%  N      - Flussvektor
%--------------------------------------------------------------------------

% ... Speicher reservieren
N = zeros(sizeZ,1);
% ... Spannungen
N(1:ntens) = - sqrt(3/2)*C*n;
% ... plastische Dehnungen
N(ntens+1:2*ntens) = sqrt(3/2)*n;
% ... Backstresstensoren
aline = sqrt( 3/2 * sum( (P2*ALPHA) .* ALPHA) );                           % Normen Backstress
aline(aline==0)=1e-40;    
hvar = aline./r; hvar(hvar>1) = 1;                                         % Hilfsvariable
hvar2 = (3/2/Y./aline) .* (B'*P2*ALPHA);                                   % Hilfsvariable
hvar2 = 0.5 * (hvar2 + abs(hvar2));                                        % Maccauley Klammer
for i = 1 : M
    N(ntens*(2+i-1)+1 : ntens*(2+i)) = c(i)*r(i)/Y * B - ...
                                   c(i)*hvar(i)^chi(i)*ALPHA(:,i)*hvar2(i);
end % Ende Schleife Backstresstensoren
end % Ende Aufstellen des Flussvektors

% -------------------------------------------------------------------------
% Aufstellen statische Erholung
% -------------------------------------------------------------------------
function M = ErholungsVektor(C,dEPS,sizeZ,ntens)
% INPUT:
% C     - Elastizitätstensor
% dEPS  - Inkrement Dehnungen
% sizeZ - Größe des Vektors
% ntens - Anzahl Tensorkomponetens 
% OUTPUT:
% M     - Erholungsvektor
%--------------------------------------------------------------------------

% ... Speicher
M = zeros(sizeZ,1);
% ... Spannungen 
M(1:ntens) = - C * dEPS;

end % Ende Aufstellen des Erholungsvektors

% -------------------------------------------------------------------------
% Aufstellen des Residuumvektors
% -------------------------------------------------------------------------
function R = Residuum(dEPS,...
                      SIG0,EPSP0,ALPHA0,...
                      SIG1,EPSP1,ALPHA1,dp1,...
                      ntens,C,P,P2,...
                      Y,M,c,r,chi)
% INPUT:
% dEPS    - Vorgegebenes Dehungsinkrement
%   SIG   - Spannungen
% EPSP    - plastische Dehnungen
% ALPHA   - Backstress
% dp1     - Inkrement der plastischen Bogenlänge
%         - _0 = Startwert           _1 = Zielwert
%  ntens  - Anzahl Tensorkomponenten
%   Y     - Radius Fließfläche
%   M     - Anzahl Backstresstensoren
% c,r,chi - Parameter OhnoWang Modell
%  C      - Elastizitätstensor
%  P      - Abbildung Spannung auf Spannungsdeviator
%  P2     - Skalarprodukt Spannungsdeviatoren & Verdopplung Von
%               Nebendiagonalelementen bei Dehnungstensoren
%
% OUTPUT:
%   dRdY - numerische Ableitung des Residuums
%
% NOTATION:
%
%       -X_n+1 + X_n + dp*N_n+1 - M_n+1
% R =   F_n+1
%
%      X - Zustandsvaroablen 
%      Y - [ X dpn+1]^T
%      N - Flussvektor
%      M - Erholungsvektor
%
%--------------------------------------------------------------------------
% Anzahl Zustandsvariablen
sizeZ = ntens*(2+M)+1;                                                     % Länge Zustandsvariablen
% ... Effektive Spannung
B1 = P * SIG1 - sum(ALPHA1,2);
% ... Normalenvektor
n1 = sqrt(3/2)*P*P2*B1./Y;
% ... Flussvektor
N1 = FlussVektor(n1,B1,ALPHA1,C,P2,Y,M,c,r,chi,sizeZ-1,ntens);
% ... Erholungsvektor
M1 = ErholungsVektor(C,dEPS,sizeZ-1,ntens);
% ... Fließfläche
F1 = MisesFF(B1,Y,P2);
% ... Residuum
dZVAR = [SIG0-SIG1;EPSP0-EPSP1;reshape(ALPHA0-ALPHA1,ntens*M,1)];
R = [ dZVAR + dp1 * N1 - M1 ; F1];

end % Ende Aufstellen Residuumsvektor

% -------------------------------------------------------------------------
% Numerisches Ableiten mit Vorwärtsdifferenzen
% -------------------------------------------------------------------------
function dRdY = VorwaertsDifferenzen(R0,dEPS,...                                                                         % Residuum
                                     SIG0,EPSP0,ALPHA0,...
                                     SIG1,EPSP1,ALPHA1,dp1,...
                                     ntens,C,P,P2,...
                                     Y,M,c,r,chi,sizeZ)
% INPUT:
% R0      - aktuelles Residuum
% dEPS    - Vorgegebenes Dehungsinkrement
%   SIG   - Spannungen
% EPSP    - plastische Dehnungen
% ALPHA   - Backstress
% dp1     - Inkrement der plastischen Bogenlänge
%         - _0 = Startwert           _1 = Zielwert
%  ntens  - Anzahl Tensorkomponenten
%   Y     - Radius Fließfläche
%   M     - Anzahl Backstresstensoren
% c,r,chi - Parameter OhnoWang Modell
%  C      - Elastizitätstensor
%  P      - Abbildung Spannung auf Spannungsdeviator
%  P2     - Skalarprodukt Spannungsdeviatoren & Verdopplung Von
%               Nebendiagonalelementen bei Dehnungstensoren
% sizeZ   - Anzahl Zustandsvariablen
% OUTPUT:
%   dRdY - numerische Ableitung des Residuums
%--------------------------------------------------------------------------

% ... Speicher reservieren
dRdY = zeros(sizeZ,sizeZ);                                                 % Speicher für die ABleitung
j = 1;                                                                     % Zeiger auf aktuelle Spalte in Ableitung
% ... Spannungen
for i = 1:ntens
    % ... bestimme Schrittweite
    vz = sign(SIG1(i)); if vz == 0, vz = 1; end
    epsilon = 1e-8 * vz * max( [abs(SIG1(i)), 1] );
    % ... Manipulierter Wert
    SIG2 = SIG1;
    SIG2(i) = SIG1(i) + epsilon;
    % ... Differenz der Residuuen
    dR = Residuum( dEPS,...                                                                         % Residuum
                   SIG0,EPSP0,ALPHA0,...
                   SIG2,EPSP1,ALPHA1,dp1,...
                   ntens,C,P,P2,...
                   Y,M,c,r,chi) - R0;
   % ... Ableitung
   dRdY(:,j) = dR./epsilon;
   % ... Inkrement Spaltenzeiger 
   j = j + 1;
end
% ... plastische Dehnungen
for i = 1:ntens
    % ... bestimme Schrittweite
    vz = sign(EPSP1(i)); if vz == 0, vz = 1; end
    epsilon = 1e-8 * vz * max( [abs(EPSP1(i)), 0.01] );
    % ... Manipulierter Wert
    EPSP2 = EPSP1;
    EPSP2(i) = EPSP2(i) + epsilon;
    % ... Differenz der Residuuen
    dR = Residuum( dEPS,...                                                                         % Residuum
                   SIG0,EPSP0,ALPHA0,...
                   SIG1,EPSP2,ALPHA1,dp1,...
                   ntens,C,P,P2,...
                   Y,M,c,r,chi) - R0;
   % ... Ableitung
   dRdY(:,j) = dR./epsilon;
   % ... Inkrement Spaltenzeiger 
   j = j + 1;
end
% ... Backstresstensoren
for k = 1 : M
    for i = 1:ntens
        % ... bestimme Schrittweite
        vz = sign(ALPHA1(i,k)); if vz == 0, vz = 1; end
        epsilon = 1e-8 * vz * max( [abs(ALPHA1(i,k)), 1] );
        % ... Manipulierter Wert
        ALPHA2 = ALPHA1;
        ALPHA2(i,k) = ALPHA2(i,k) + epsilon;
        % ... Differenz der Residuuen
        dR = Residuum( dEPS,...                                                                         % Residuum
                       SIG0,EPSP0,ALPHA0,...
                       SIG1,EPSP1,ALPHA2,dp1,...
                       ntens,C,P,P2,...
                       Y,M,c,r,chi) - R0;
        % ... Ableitung
        dRdY(:,j) = dR./epsilon;
        % ... Inkrement Spaltenzeiger 
        j = j + 1;
    end % Ende Schleife Tensorkomponenten der Backstresstensoren
end % Ende Schleife Backstresstensoren
% ... Inkrement plastische Bogenlänge
% ... bestimme Schrittweite
epsilon = 1e-8 * max( [abs(dp1), 1e-5] );
% ... Manipulierter Wert
dp2 = dp1 + epsilon;
% ... Differenz der Residuuen
dR = Residuum( dEPS,...                                                                         % Residuum
    SIG0,EPSP0,ALPHA0,...
    SIG1,EPSP1,ALPHA1,dp2,...
    ntens,C,P,P2,...
    Y,M,c,r,chi) - R0;
% ... Ableitung
dRdY(:,j) = dR./epsilon;
end % Ende Vorwärtsdifferenzen



