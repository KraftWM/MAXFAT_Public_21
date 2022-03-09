function [ZVARneu,DEP,CEP] = OWT(ntens, ndi, ink, ZVAR, ink_flag, parameter)
%
% (O)hno (W)ang (T)anaka
%
% kinematische Verfestigung:
% Aus Ohno et al. 1993 KINEMATIC HARDENING RULES WITH CRITICAL
% STATE OF DYNAMIC RECOVERY, PART I
%
% nichtproportionale Verfestigung:
% Vorbild aus Döring mit Tanaka Parameter
% 
%
% Implentierung nach Döring
%
%   INPUT:
%         ntens     -> Anzahl Tensorkomponten (in Spannungen)
%         ndi       -> Anzahl Diagonalelemente, mehr kkommentare
%         ink       -> Belastungsinkrement
%         ZVAR      -> Zustandsvariablen [eps;epsp;alphai;p,Q,beta,q,A,C] bein spansteu
%                                        [sig;epsp;alphai;p,Q,beta,q,A,C] bein dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         parameter -> Materialparameter des Modells
%                      [E,nu,...
%                       c_i,r_i,chi_i,...
%                       gamma,gamma_np,gamma_a,gamma_c,...
%                       Qnpmax,...
%                       eta,omega,cg,...
%                       r0]
%
%
%         ZVAR(1:ntens)   = eps                                                   -> Gesamtdehnung
%         ZVAR(1:ntens)   = sig                                                   -> Spannungen
%         ZVAR(ntens+1:2*ntens)  = epsp                                           -> plastische Dehnungen
%         ZVAR(2*ntens+1:(2+M)*ntens) = alpha_i                                   -> Teilbackstresstensoren
%         ZVAR((2+M)*ntens+1) =  p                                                -> plastische Bogenlänge, hier dp = sqrt(dep:dep)
%         ZVAR((2+M)*ntens+2) =  Q                                                -> Zusätzlicher Radius FF für NPV
%         ZVAR((2+M)*ntens+3:(3+M)*ntens+2) = beta                                -> (Backstraintensor) Mittelpunkt Gedächtnisfläche im Dehnungsraum
%         ZVAR((3+M)*ntens+3) = q                                                 -> Radius Dehnungsgedächtnisfläche
%         ZVAR((3+M)*ntens+4) = A                                                 -> Nichtproportionaalitätskennwert
%         ZVAR((3+M)*ntens+5 : (3+M+(ntens+1)/2)*ntens+4) = CT                    -> Nichtproportionalitätstensor nach Tanaka
%
%
%    OUTPUT:
%        ZVARneu -> neue zustandsvariablen nach Lastinkrement
%    
%    NOTATIONEN:
%    Spannungen                 Verzerrungen
%
%           sig_11                     eps_11
%    sig =  sig_22               eps = eps_22
%           sig_12                     2eps_12
%
%    Tanaka Tensor
%     -> Symmetrischer Tensor 4. Stufe (Dehnungscharakter)
%     -> 21 unabhängige Komponenten (bei ntens = 6) (num Komponenten:
%        ntens*(ntens+1)/2
%     -> Speicherrichtung entlang der Diagonalen
%            Tensorkomp. in Voigt Notation                      Für Zustandsvariablen
%        ( 1111   1122   1133  2*1112  2*1113  2*1123 )       ( 1  7  12  16  19  21)
%        (        2222   2233  2*2212  2*2213  2*2223 )       (    2   8  13  17  20)
%   CT = (               3333  2*3312  2*3312  2*3323 )  =    (        3   9  14  18)
%        (                     4*1212  4*1213  4*1223 )       (            4  10  15)
%        (                             4*1313  4*1323 )       (                5  11)
%        (                                     4*2323 )       (                    6)
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: Juni 2021                                                    |
%  ----------------------------------------------------------------------

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------

% Anzahl Backstress
M = (length(parameter)-11)/3; % TODO
% Elastizitätskonstanten
E = parameter(1);
nu = parameter(2);
r0 = parameter(end); % TODO

%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% Steifigkeit
C = elast_steifigkeit(E,nu,ntens,ndi);
% statische Matrizen
[P, P_line] = set_maps(ntens,ndi);

%--------------------------------------------------------------------------
%                  Zustandsvariablen auslesen                                 
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
if ink_flag == 0                 % spansteu
    eps = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
    sig = C * (eps - epsp);
    dsig_tr = ink;               % elastischer trial step
elseif ink_flag == 1             % dehnsteu
    sig = ZVAR(1:ntens);
    dsig_tr = C * ink;           % elastischer trial step
else
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end
% Backstresstensoren
alpha = reshape( ZVAR(2*ntens+1:(M+2)*ntens) , ntens, M);
% Radius FF
Q = ZVAR((2+M)*ntens+2);
r = r0 + Q;


%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------

% Trial Inkrement der Spannungen 
if ntens == 1
    msg = ['Modell nicht implementiert. Im einachsigen gibts keine ',...
           'nichtproportionale Verfestigung. Nimm doch einfach das ',...
           'Ohno Wang Modell Arschloch !!!'];
    error(msg)
elseif ntens == 2 % sigma-tau
    s = P .* sig;                                                          % Span dev
    ds = P .* dsig_tr;                                                     % dev des Spannungsinkrements
    s_tr = s + ds;                                                         % Trial Spannungsdeviator
    a = sum(alpha,2);                                                      % Backstress;
    beta = s_tr - a;                                                       % Relativ Spannung
    F_tr = beta' * (P_line .* beta) - 2/3 * r^2;
else
    s = P * sig;                                                           % Span dev
    ds = P * dsig_tr;                                                      % dev des Spannungsinkrements
    s_tr = s + ds;                                                         % Trial Spannungsdeviator
    a = sum(alpha,2);                                                      % Backstress;
    beta = s_tr - a;                                                       % Relativ Spannung
    F_tr = beta' * P_line * beta - 2/3 * r^2;                              % Trial Fließfunktion
end
FTOL = 1e-7;                                                               % Toleranz für Abweichungen F ~= 0


%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------

if F_tr < FTOL   % Trial step völlig elastische
    
    if ink_flag == 0 % spansteu;
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        eps = eps + D * ink;
        ZVARneu = [eps;ZVAR(ntens+1:end)];
    elseif ink_flag == 1 % dehnsteu
        sig = sig + dsig_tr;
        ZVARneu = [sig;ZVAR(ntens+1:end)];
    end
    
    % Falls tangentiale Steifigkeit gebraucht wird
    if nargout == 2 % Berechne elastisch plastische nachgiebigkeit
%         CEP = elast_steifigkeit(E,nu,ntens,ndi);
          DEP = elast_nachgiebigkeit(E,nu,ntens,ndi);
    elseif nargout == 3
          CEP = elast_steifigkeit(E,nu,ntens,ndi);
          DEP = elast_nachgiebigkeit(E,nu,ntens,ndi);
    end


%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%-------------------------------------------------------------------------- 

else % Ermittle elastischen Anteil xel
    
    % 1. ermittle elastischen Anteil
%     xel = elastink(s,a,r0,ds,ndi);
    xel = elastink2(s,a,P_line,r,ds,FTOL);
    
    % 2. ausführen elastischen Teil
    if ink_flag == 0 % spansteu;
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        eps = eps + D * (xel*ink);
        ZVAR = [eps;ZVAR(ntens+1:end)];
    elseif ink_flag == 1 % dehnsteu
        sig = sig + xel*dsig_tr;
        ZVAR = [sig;ZVAR(ntens+1:end)];
    end
    
    % 3. plastischer Teil des inkrements
    ink = (1-xel)*ink;
    
    % Integration je nach Spannungszustand
    options = [];
    if ntens == 6 % 3D
        
        % nachg.
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % abbildungen
        [P, P_hat, P_tilde] = set_maps(ntens,ndi);
        
        % Integration
        
%         % Matlab Version
%         [~,ZVARneu] = rk87(@materialmodell3d,[0,1], ZVAR, options, ink,...
%                         ink_flag, parameter, C, D, P, P_hat, P_tilde);

        % Mex Version
        [~,ZVARneu] = rk87(@CoderVersion_OWT_3D_mex,[0,1], ZVAR, options,...
                        M, ink, ink_flag, parameter, C, D, P, P_hat, P_tilde);
                    
    elseif ntens == 3 && ndi == 2 % ESZ
        
        % nachgiebigkeitsmatrix
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % statische Matrizen
        [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
        
        % Integration
        
        % matlab Version
%         [~,ZVARneu] = rk87(@materialmodellESZ,[0,1], ZVAR, options, ink,...
%                         ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);

        % Mex Version
        [~,ZVARneu] = rk87(@CoderVersion_OWT_ESZ_mex,[0,1], ZVAR, options,...
                        M, ink, ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);
    
    elseif ntens == 2 && ndi == 1 % sigma-tau
        
        % nachgiebigkeitsmatrix
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % statische Matrizen
        [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
        
        % Integration
        
        % matlab Version
%         [~,ZVARneu] = rk87(@materialmodellST,[0,1], ZVAR, options, ink,...
%                         ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);
        
%          Mex Version
        [~,ZVARneu] = rk87(@CoderVersion_OWT_ST_mex,[0,1], ZVAR, options,...
                        M, ink, ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);
                
        
    elseif ntens == 1 % 1D
        
        msg = ['Modell nicht implementiert. Im einachsigen gibts keine ',...
           'nichtproportionale Verfestigung. Nimm doch einfach das ',...
           'Ohno Wang Modell Arschloch !!!'];
        error(msg)
                    
    end
    
    % nur letzten Schritt ausgeben
    ZVARneu = ZVARneu(end,:)'; 
    
    % Testen der Konsistenzbedingung
    if ntens~=1
        ZVARneu = konsistenzbedingung(ZVARneu,M,P,P_line,ntens,r0,C,D,ink_flag);
    end
    
    % Falls tangentiale Steifigkeit gebraucht wird
    if nargout == 2 % berechne Elastisch plastische steifigkeitsmatrix
        msg = 'noch nicht implementiert';
        error(msg);
%         C = elast_steifigkeit(E,nu,ntens,ndi);
%         DEP = tangnach_ohnowang(ZVARneu,parameter,ink_flag,ntens,ndi,C,D);
    elseif nargout == 3
        msg = 'noch nicht implementiert';
        error(msg);
%         C = elast_steifigkeit(E,nu,ntens,ndi);
%         DEP = tangnach_ohnowang(ZVARneu,parameter,ink_flag,ntens,ndi,C,D);
%         CEP = tangsteif_ohnowang(ZVARneu,parameter,ink_flag,ntens,ndi,C);
    end
    
end % Ende Unterscheidung Trial Step
end % Ende Hauptfunktion



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Prüfe die konsistenzbedingung                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = konsistenzbedingung(X,M,P,P_line,ntens,r,C,D,ink_flag)
% Window Ausgabe
% fprintf('Zustandsvariablen Original:\n')
% fprintf('%.32d \n',X)

% Radius FF
r = r + X((2+M)*ntens+2);

% auslesen zustände
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
    % epsp = X(ntens+1:2*ntens);
end
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
a = sum(alpha,2);


% Relativspannung
% Fließfläche
if ntens == 1
    s = P * sig;
    beta = s - a;
    SV = beta^2;
elseif ntens == 2
    s = P .* sig;
    beta = s - a;
    SV = 1.5 * beta' * (P_line .* beta);
else
    s = P * sig;
    beta = s - a;
    SV = 1.5 * beta' * P_line * beta;
end
% F1 = SV - r^2;
% fprintf('Überspannung Orig: %.4d    ',F1)

% Prüfe Konsistenzbedingung
if SV - r^2 > 1e-11
    
    % Fehler/Warnung wenn zu stark verletzt
    if SV - r^2 > 1e-1
        msg = 'Konsistenzbedinung verletzt';
        error(msg)
    elseif SV - r^2 > 1e-4
        msg = 'Konsistenzbedinung geringfügig verletzt';
        warning(msg)
    end

    % Korrigiere Spannungsdeviator (Rel. Span. radial auf FF zurückprojezieren)
    fak = r/sqrt(SV);
    if ntens == 6
        % ... 3D nur Deviator Korrigieren
        s = a + fak*beta;
        sh = (sig(1)+sig(2)+sig(3))/3;
%         eins = [1; 1; 1; 0; 0; 0];
        % Korrigierte Spannung 
%         sig = s + [sh; sh; sh; 0; 0; 0];
        sig = s;
        sig(1) = sig(1) + sh;
        sig(2) = sig(2) + sh;
        sig(3) = sig(3) + sh;
    elseif ntens == 3
        % ... ESZ nur Deviator korriegieren 
        s = a + fak*beta;
        % korrigiere Spannungen
        sig = [2 1 0; 1 2 0; 0 0 1] * s;
%         sh = (sig(1)+sig(2))/3;
%         eins = [1; 1; 0];
    elseif ntens == 2
        % ... sig-tau nur Deviator korrigieren
        s = a + fak*beta;
        % korrigiere Spannungen
        sig = [1.5*s(1);s(2)];
    elseif ntens == 1
        % ... 1D Spannung korrigieren
        sh = 0;
        s = a + fak*beta;
        eins = 0;
        % Korrigierte Spannung 
        sig = sh*eins + s;
    end
    
    

    % Test zum debuggen
%     beta = s - a;
%     SV = 1.5 * beta' * P_line * beta;
%     F2 = SV - r^2;
%     fprintf('Überspannung Korr1: %.4d \n',F2)
%     fprintf('Spanndev korr:\n')
%     fprintf('%.32d\n',s)

    % Korrigierte Zustandsgrößen
    if ink_flag == 0 % Spansteu.
        epse = D*sig;
        eps = epse+epsp;
        X(1:ntens) = eps;
%         sig = C * (eps-epsp);
%         s = P * sig ;
%         fprintf('Spanndev korr3:\n')
%         fprintf('%.32d\n',s)
%         fprintf('\n\n')
%         beta = P*sig - a;
%         SV = 1.5 * beta' * P_line * beta;
%         F2 = SV - r^2;
%         fprintf('Überspannung Korr2: %.4d \n\n',F2)
    elseif ink_flag == 1 % Dehnsteu
        X(1:ntens) = sig;
%         if ntens == 2
%             beta = P.*sig - a;
%             SV = 1.5 * beta' * (P_line .* beta);
%         else
%             beta = P*sig - a;
%             SV = 1.5 * beta' * P_line * beta;
%         end
%         F2 = SV - r^2;
%         fprintf('Überspannung Korr2: %.4d \n\n',F2)
    end

    % fprintf('Zustandsvariablen Korrigiert:\n')
    % fprintf('%.32d \n',X)
end
% fprintf('\n')
end % Ende Prüfen Konsitenzbedingung





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 3D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dZVAR = ...
          materialmodell3d(~, ZVAR, ink, ink_flag, parameter , C, D, ...
                              P, P_hat, P_tilde)
% Konkretes Materialmodell für 3D Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
% Materialmodell:
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_hat ...    -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------


% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 6;                                                                 % Tensorkomponenten
M = (length(parameter)-11)/3;   % TODO                                            % Anzahl TBST

% elastische 

r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);
% NP Verfestigung
gamma = parameter(3+3*M);
gamma_np = parameter(4+3*M);
gamma_a = parameter(5+3*M);
gamma_c = parameter(6+3*M);
Qnpmax = parameter(7+3*M);
eta = parameter(8+3*M);
omega = parameter(9+3*M);
cg = parameter(10+3*M);

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
delta = 1e-40;

% -------------------------------------------------------------------------
%                        Zustandsvariablen
% -------------------------------------------------------------------------

% auslesen zustände
if ink_flag == 0 % Spansteu.
    eps = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
end
alpha = reshape(ZVAR(2*ntens+1:(M+2)*ntens),ntens,M);
% p = ZVAR((M+2)*ntens+1);
Q = ZVAR((M+2)*ntens+2);
beta = ZVAR((2+M)*ntens+3:(3+M)*ntens+2);
q = ZVAR((3+M)*ntens+3);
A = ZVAR((3+M)*ntens+4);
CT = ZVAR((3+M)*ntens+5 : (3+M+(ntens+1)/2)*ntens+4);

% -------------------------------------------------------------------------
%                        Radius der Fließfläche
% -------------------------------------------------------------------------
r = r0 + Q;


% -------------------------------------------------------------------------
%                   Normale an die Fließfläche
% -------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r).*P_hat*(s-a);
transn = n.';

% -------------------------------------------------------------------------
%                   Ableitung der Teilbacksztresstensoren
% -------------------------------------------------------------------------

% normen der Teilbackstresstensoren 
norm_ai = sqrt( sum( (P_hat*alpha) .* alpha) );
norm_ai(norm_ai == 0) = delta;

% Hilfsvariablen 
Li = alpha./norm_ai;
var1 = c_i .* r_i;
var2 = P_tilde * n;
var3 = norm_ai./r_i;
var3(var3>1) = 1;
var4 = transn * Li;
var4 = 0.5 * (var4 + abs(var4));
% Schleife über alle Backstresstensoren (aus irgend einem Grund ist
% Schleife schneller als vektorisiert)
% init ableitungen
dalpha_dp = zeros(ntens,M);
for ii = 1 : M
    
    % backstress
    dalpha_dp(:,ii) = var1(ii) * (var2 - (var3(ii)).^(chi_i(ii)+1).*var4(ii).*Li(:,ii));
    
end
% Ableitungen der gesamtbackstresstensoren
da_dp = sum(dalpha_dp,2);

% -------------------------------------------------------------------------
%                   Gedächtnissfläche
% -------------------------------------------------------------------------
% Effektive plastische Dehnung
effstrain = epsp-beta;
norm2es = effstrain'*P_tilde*effstrain;

% Gedächtnissfläche
g = norm2es - q^2;

% Hilfsvariable ( H(g) )
Hg = 0.5 * (sign(g) + 1);

% Normale an die gedächtnissfläche
if norm2es == 0
    norm2es = delta;
end
if Hg == 0 
    % ... ||n*|| < 1
%     ng = effstrain./sqrt(norm2es);
    ng = effstrain./q;%effstrain./sqrt(norm2es);
    dqdp = - (cg * q)^omega; 
else
    % ... ||n*|| = 1
    if norm2es > delta
        ng = effstrain./sqrt(norm2es);
    else
        ng = n;
    end
    dqdp = eta;
end


% Hilfsvariabel <ne:n>
nn = transn * P_tilde * ng;
nn = 0.5 * (nn + abs(nn));

% Isotrope Verfestigung
dqdp = dqdp*nn;

% Kinematische Verfestigung
dbetadp = (1-eta) * Hg * nn * ng;

% -------------------------------------------------------------------------
%                   Nichtproportionale Verfestigung
% -------------------------------------------------------------------------
Qnpinf = Qnpmax * (1 - exp(- gamma_np * q));
Qnp = A*Qnpinf;
dQdp = gamma * ( Qnp - Q );

% -------------------------------------------------------------------------
%                   plastischer Tangentenmodul
% -------------------------------------------------------------------------

h = transn * da_dp + sqrt(2/3)* dQdp;

% -------------------------------------------------------------------------
%                   Inkrement plastische Bogenlänge
% -------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (1/h) * (transn * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = ( transn * (C * ink)) / ...
         (h +  transn * ( C * n));
else % Fehler
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end

% -------------------------------------------------------------------------
%                   Inkremente Zustandsvariableb
% -------------------------------------------------------------------------

depsp= dp .* n;
dalpha= dalpha_dp .* dp;
dQ = dQdp * dp;
dq = dqdp * dp;
dbeta = dbetadp * dp;
[dA,dCT] = inkTanaka(ntens,A,CT,n,gamma_a,gamma_c,dp);

% -------------------------------------------------------------------------
%                   zusammenfassen der Inkremente
% -------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
dZVAR = [dout; depsp ; reshape(dalpha,ntens*M,1); dp; dQ; dbeta; dq; dA; dCT];   

end % Ende Modell 3D



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen Ebener Spannungszustand             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dZVAR = ...
          materialmodellESZ(~, ZVAR, ink, ink_flag, parameter, C, D, P, ...
                            P_line, P_hat, A, P_check)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(parameter)-11)/3;   % TODO                                            % Anzahl TBST

% elastische 

r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);
% NP Verfestigung
gamma = parameter(3+3*M);
gamma_np = parameter(4+3*M);
gamma_a = parameter(5+3*M);
gamma_c = parameter(6+3*M);
Qnpmax = parameter(7+3*M);
eta = parameter(8+3*M);
omega = parameter(9+3*M);
cg = parameter(10+3*M);

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
delta = 1e-40;                                                             % numerisch 0

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% auslesen zustände
if ink_flag == 0 % Spansteu.
    eps = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
end
alpha = reshape(ZVAR(2*ntens+1:(M+2)*ntens),ntens,M);
% p = ZVAR((M+2)*ntens+1);
Q = ZVAR((M+2)*ntens+2);
beta = ZVAR((2+M)*ntens+3:(3+M)*ntens+2);
q = ZVAR((3+M)*ntens+3);
ANP = ZVAR((3+M)*ntens+4);
CT = ZVAR((3+M)*ntens+5 : (3+M+(ntens+1)/2)*ntens+4);


% -------------------------------------------------------------------------
%                        Radius der Fließfläche
% -------------------------------------------------------------------------
r = r0 + Q;


%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r).* P_hat * (s-a);
transn = n.';


%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Normen der Teilbackstresstensoren
norm_ai = sqrt( sum( (P_line*alpha) .* alpha ) );
norm_ai(norm_ai == 0) = delta;                                             % Abfangen von 0

% Hilfsvariablen 
Li = alpha./norm_ai;
var1 = c_i .* r_i;
var2 = A * n;
var3 = norm_ai./r_i;
var3(var3>1) = 1;
var4 = transn * P_check * Li;
var4 = 0.5 * (var4 + abs(var4));
% Schleife über alle Backstresstensoren (aus irgend einem Grund ist
% Schleife schneller als vektorisiert)
% init ableitungen
dalpha_dp = zeros(ntens,M);
for ii = 1 : M
    
    % backstress
    dalpha_dp(:,ii) = var1(ii) * (var2 - ...
           (var3(ii)).^(chi_i(ii)+1).*var4(ii).*Li(:,ii));
    
    % backstress
%     dalpha_dp(:,ii) = h_i(ii) * (var2 - ...
%                              (var3(ii)).^(chi_i(ii)+1).*var4(ii).*Li(:,ii));
    
end

% Ableitung gesamt backstress
da_dp = sum(dalpha_dp,2);

% -------------------------------------------------------------------------
%                   Gedächtnissfläche
% -------------------------------------------------------------------------
% Effektive plastische Dehnung
effstrain = epsp-beta;
norm2es = effstrain'*(A*P_check)*effstrain;

% Gedächtnissfläche
g = norm2es - q^2;

% Hilfsvariable ( H(g) )
Hg = 0.5 * (sign(g) + 1);

% Normale an die gedächtnissfläche
if norm2es == 0
    norm2es = delta;
end
if Hg == 0 
    % ... ||n*|| < 1
%     ng = effstrain./sqrt(norm2es);
    ng = effstrain./q;%effstrain./sqrt(norm2es);
    dqdp = - (cg * q)^omega; 
else
    % ... ||n*|| = 1
    if norm2es > delta
        ng = effstrain./sqrt(norm2es);
    else
        ng = n;
    end
    dqdp = eta;
end


% Hilfsvariabel <ne:n>
nn = transn * (A*P_check) * ng;
nn = 0.5 * (nn + abs(nn));

% Isotrope Verfestigung
dqdp = dqdp*nn;

% Kinematische Verfestigung
dbetadp = (1-eta) * Hg * nn * ng;

% -------------------------------------------------------------------------
%                   Nichtproportionale Verfestigung
% -------------------------------------------------------------------------
Qnpinf = Qnpmax * (1 - exp(- gamma_np * q));
Qnp = ANP*Qnpinf;
dQdp = gamma * ( Qnp - Q );



%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = transn * P_check * da_dp + sqrt(2/3) * dQdp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (1/h) * (transn * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = (transn * (C * ink)) / ...
         (h + transn * ( C * n));
else % Fehler
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp= dp .* n;
dalpha= dalpha_dp .* dp;
dQ = dQdp * dp;
dq = dqdp * dp;
dbeta = dbetadp * dp;
[dA,dCT] = inkTanaka(ntens,ANP,CT,n,gamma_a,gamma_c,dp);

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end

dZVAR = [dout; depsp ; reshape(dalpha,ntens*M,1); dp; dQ; dbeta; dq; dA; dCT];   

end % Ende Modell ESZ



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen sigma-tau Spannungszustand          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dZVAR = ...
          materialmodellST(~, ZVAR, ink, ink_flag, parameter, C, D, P, ...
                            P_line, P_hat, A, P_check)
% Konkretes Materialmodell für sigma-tau Spannungszustände, gibt bei 
% vorgabe eines Lastinkrementes die Inkremente der inneren Variablen 
% zurück. Dabei wird % angenommen, das jedes übergebene Inkrement 
% elastisch-plastische Deformationen hervorruft.
%
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2021 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

% Zuweisen der Materialparameter
ntens = 2;                                                                 % Tensorkomponenten
M = (length(parameter)-11)/3;                                              % Anzahl TBST

% elastische 
r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);
% NP Verfestigung
gamma = parameter(3+3*M);
gamma_np = parameter(4+3*M);
gamma_a = parameter(5+3*M);
gamma_c = parameter(6+3*M);
Qnpmax = parameter(7+3*M);
eta = parameter(8+3*M);
omega = parameter(9+3*M);
cg = parameter(10+3*M);

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
delta = 1e-40;

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
if ink_flag == 0 % Spansteu.
    eps = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
end
% backstresstensoren
alpha = reshape(ZVAR(2*ntens+1:(M+2)*ntens),ntens,M);
% p = ZVAR((M+2)*ntens+1);
Q = ZVAR((M+2)*ntens+2);
beta = ZVAR((2+M)*ntens+3:(3+M)*ntens+2);
q = ZVAR((3+M)*ntens+3);
ANP = ZVAR((3+M)*ntens+4);
CT = ZVAR((3+M)*ntens+5 : (3+M+(ntens+1)/2)*ntens+4);

% -------------------------------------------------------------------------
%                        Radius der Fließfläche
% -------------------------------------------------------------------------
r = r0 + Q;

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P .* sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r).* P_hat .* (s-a);
transn = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Normen der Teilbackstresstensoren
norm_ai = sqrt( sum( (P_line.*alpha) .* alpha ) );
norm_ai(norm_ai == 0) = delta;                                             % Abfangen von 0
% % Hilfsgrößen
Li = alpha./norm_ai;
var1 = c_i .* r_i;
var2 = A .* n;
var3 = norm_ai./r_i;
var3(var3>1) = 1;
var4 = transn * (P_check .* Li);
var4 = 0.5 * (var4 + abs(var4));
% Schleife über alle Backstresstensoren (aus irgend einem Grund ist
% Schleife schneller als vektorisiert)
% init ableitungen
dalpha_dp = zeros(ntens,M);
for ii = 1 : M
    
    % backstress
    dalpha_dp(:,ii) = var1(ii) * (var2 - ...
           (var3(ii)).^(chi_i(ii)+1).*var4(ii).*Li(:,ii));

end

% Ableitung gesamt backstress
da_dp = sum(dalpha_dp,2);

% -------------------------------------------------------------------------
%                   Gedächtnissfläche
% -------------------------------------------------------------------------
% Effektive plastische Dehnung
effstrain = epsp-beta;
norm2es = effstrain' * ((A.*P_check) .* effstrain);

% Gedächtnissfläche
g = norm2es - q^2;

% Hilfsvariable ( H(g) )
Hg = 0.5 * (sign(g) + 1);

% Normale an die gedächtnissfläche
if norm2es == 0
    norm2es = delta;
end
if Hg == 0 
    % ... ||n*|| < 1
%     ng = effstrain./sqrt(norm2es);
    ng = effstrain./q;%effstrain./sqrt(norm2es);
    dqdp = - (cg * q)^omega; 
else
    % ... ||n*|| = 1
    if norm2es > delta
        ng = effstrain./sqrt(norm2es);
    else
        ng = n;
    end
    dqdp = eta;
end


% Hilfsvariabel <ne:n>
nn = transn * ((A.*P_check) .* ng);
nn = 0.5 * (nn + abs(nn));

% Isotrope Verfestigung
dqdp = dqdp*nn;

% Kinematische Verfestigung
dbetadp = (1-eta) * Hg * nn * ng;

% -------------------------------------------------------------------------
%                   Nichtproportionale Verfestigung
% -------------------------------------------------------------------------
Qnpinf = Qnpmax * (1 - exp(- gamma_np * q));
Qnp = ANP*Qnpinf;
dQdp = gamma * ( Qnp - Q );

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = transn * (P_check .* da_dp) + sqrt(2/3) * dQdp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (1/h) * (transn * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = (transn * (C * ink)) / ...
         (h + transn * ( C * n));
else % Fehler
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp= dp .* n;
dalpha= dalpha_dp .* dp;
dQ = dQdp * dp;
dq = dqdp * dp;
dbeta = dbetadp * dp;
[dANP,dCT] = inkTanaka(ntens,ANP,CT,n,gamma_a,gamma_c,dp);

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
dZVAR = [dout; depsp ; reshape(dalpha,ntens*M,1); dp; dQ; dbeta; dq; dANP; dCT];   

end % Ende Modell Sigma-tau Spannungszustand