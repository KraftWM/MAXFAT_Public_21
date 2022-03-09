function [X_neu,DEP,CEP] = ohnowang(ntens, ndi, ink, X, ink_flag, parameter)
% Materialmodell:
% Aus Ohno et al. 1993 KINEMATIC HARDENING RULES WITH CRITICAL
% STATE OF DYNAMIC RECOVERY, PART I
%
% Implentierung nach Döring
%
%   INPUT:
%         ntens     -> Anzahl Tensorkomponten (in Spannungen)
%         ndi       -> Anzahl Diagonalelemente, mehr kkommentare
%         ink       -> Belastungsinkrement
%         X         -> Zustandsvariablen [eps;epsp;alphai;p] bein spansteu
%                                        [sig;epsp;alphai;p] bein dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         parameter -> Materialparameter des Modells
%                      [E,nu,c_i,r_i,chi_i,r0]
%
%
%    OUTPUT:
%        X_neu -> neue zustandsvariablen nach Lastinkrement
%    
%    NOTATIONEN:
%    Spannungen                 Verzerrungen
%
%           sig_11                     eps_11
%    sig =  sig_22               eps = eps_22
%           sig_12                     2eps_12
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: januar 2019                                                  |
%  ----------------------------------------------------------------------

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------

% Anzahl Backstress
M = (length(parameter)-3)/3;
% Elastizitätskonstanten
E = parameter(1);
nu = parameter(2);
r0 = parameter(end);

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
if ink_flag == 0 % spansteu
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
    dsig_tr = ink;    % elastischer trial step
elseif ink_flag == 1 % dehnsteu
    sig = X(1:ntens);
    dsig_tr = C * ink; % elastischer trial step
else
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end
% Backstresstensoren
alpha = reshape( X(2*ntens+1:(M+2)*ntens) , ntens, M);


%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------

% Trial Inkrement der Spannungen 
if ntens == 1
    s = P * sig;                % Span dev
    ds = P * dsig_tr;           % dev des Spannungsinkrements
    s_tr = s + ds;              % Trial Spannungsdeviator
    a = sum(alpha,2);           % Backstress;
    beta = s_tr - a;                    % Relativ Spannung
    F_tr = abs(beta) - r0;
elseif ntens == 2 % sigma-tau
    s = P .* sig;                                                           % Span dev
    ds = P .* dsig_tr;                                                      % dev des Spannungsinkrements
    s_tr = s + ds;                                                         % Trial Spannungsdeviator
    a = sum(alpha,2);                                                      % Backstress;
    beta = s_tr - a;                                                       % Relativ Spannung
    F_tr = beta' * (P_line .* beta) - 2/3 * r0^2;
else
    s = P * sig;                % Span dev
    ds = P * dsig_tr;           % dev des Spannungsinkrements
    s_tr = s + ds;              % Trial Spannungsdeviator
    a = sum(alpha,2);           % Backstress;
    beta = s_tr - a;                    % Relativ Spannung
    F_tr = beta' * P_line * beta - 2/3 * r0^2; % Trial Fließfunktion
end
FTOL = 1e-7;                        % Toleranz für Abweichungen F ~= 0


%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------

if F_tr < FTOL   % Trial step völlig elastische
    
    if ink_flag == 0 % spansteu;
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        eps = eps + D * ink;
        X_neu = [eps;X(ntens+1:end)];
    elseif ink_flag == 1 % dehnsteu
        sig = sig + dsig_tr;
        X_neu = [sig;X(ntens+1:end)];
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
    xel = elastink2(s,a,P_line,r0,ds,FTOL);
    
    % 2. ausführen elastischen Teil
    if ink_flag == 0 % spansteu;
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        eps = eps + D * (xel*ink);
        X = [eps;X(ntens+1:end)];
    elseif ink_flag == 1 % dehnsteu
        sig = sig + xel*dsig_tr;
        X = [sig;X(ntens+1:end)];
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
        
        % Matlab Version
%         [~,X_neu] = rk87(@materialmodell3d,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D, P, P_hat, P_tilde);

        % Mex Version
        [~,X_neu] = rk87(@CoderVersion_3D_OhnoWang_mex,[0,1], X, options,...
                        M, ink, ink_flag, parameter, C, D, P, P_hat, P_tilde);
                    
    elseif ntens == 3 && ndi == 2 % ESZ
        
        % nachgiebigkeitsmatrix
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % statische Matrizen
        [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
        
        % Integration
        
        % matlab Version
%         [~,X_neu] = rk87(@materialmodellESZ,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);

        % Mex Version
        [~,X_neu] = rk87(@CoderVersion_ESZ_OhnoWang_mex,[0,1], X, options,...
                        M, ink, ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);
    
    elseif ntens == 2 && ndi == 1 % sigma-tau
        
        % nachgiebigkeitsmatrix
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % statische Matrizen
        [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
        
        % Integration
        
        % matlab Version
%         [~,X_neu] = rk87(@materialmodellST,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);
        
        % Mex Version
        [~,X_neu] = rk87(@CoderVersion_OhnoWang_ST_mex,[0,1], X, options,...
                        M, ink, ink_flag, parameter, C, D,P, P_line, P_hat, A, P_check);
    
                    
    elseif ntens == 1 % 1D
        
        D = 1/C;
        
        [~,X_neu] = rk87(@materialmodell1d,[0,1], X, options, ink,...
                        ink_flag, parameter);
                    
    end
    
    % nur letzten Schritt ausgeben
    X_neu = X_neu(end,:)'; 
    
    % Testen der Konsistenzbedingung
    X_neu = konsistenzbedingung(X_neu,M,P,P_line,ntens,r0,C,D,ink_flag);
    
    
    % Falls tangentiale Steifigkeit gebraucht wird
    if nargout == 2 % berechne Elastisch plastische steifigkeitsmatrix
%         CEP = tangsteifigkeit_ohnowang(X_neu,parameter,ink_flag,ntens,ndi);
        C = elast_steifigkeit(E,nu,ntens,ndi);
        DEP = tangnach_ohnowang(X_neu,parameter,ink_flag,ntens,ndi,C,D);
    elseif nargout == 3
        C = elast_steifigkeit(E,nu,ntens,ndi);
        DEP = tangnach_ohnowang(X_neu,parameter,ink_flag,ntens,ndi,C,D);
        CEP = tangsteif_ohnowang(X_neu,parameter,ink_flag,ntens,ndi,C);
    end
    
end
end % Ende Hauptfunktion







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Prüfe die konsistenzbedingung                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = konsistenzbedingung(X,M,P,P_line,ntens,r,C,D,ink_flag)
% Window Ausgabe
% fprintf('Zustandsvariablen Original:\n')
% fprintf('%.32d \n',X)

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

% fprintf('Spanndev orig:\n')
% fprintf('%.32d\n',s)
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
% fprintf('Überspannung Orig: %.4d \n',F1)

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
%         beta = P*sig - a;
%         SV = 1.5 * beta' * P_line * beta;
%         F2 = SV - r^2;
%         fprintf('Überspannung Korr2: %.4d \n\n',F2)
    end

    % fprintf('Zustandsvariablen Korrigiert:\n')
    % fprintf('%.32d \n',X)
end
end % Ende Prüfen Konsitenzbedingung














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Elastisch Plastische tangentiale Nachgiebigkeit            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DEP = tangnach_ohnowang(X,parameter,ink_flag,ntens,ndi,CEL,DEL)
% Funktion berechnet nachgiebigkeit im elastisch-plastischen

% Abbildungen
if ntens == 6
    [P,Paa] = set_maps(ntens,ndi);
    Pna = diag(ones(1,6));
elseif ntens == 3
    [P,Paa,~,~,Pna] = set_maps(ntens,ndi);
else
    [P,Paa,Pna] = set_maps(ntens,ndi);
end
% Parameter
M = (length(parameter)-3)/3;                                               % Anzahl TBST
r0 = parameter(end);                                                       % startradius fliessfläche
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);
% Zustandsvariablen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = CEL * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
    % epsp = X(ntens+1:2*ntens);
end
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% Normale an die Fließfläche
n = sqrt(3/2) * diag([1,1,2])*(s - a) / r0;
% tangentenmodul
h = tangmod(c_i,r_i,chi_i,alpha,n,Pna,Paa);
% Tangetiale Nachgiebigkeit
DEP = DEL + (n * n')/h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Elastisch Plastische tangentiale Steifigkeit               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CEP = tangsteif_ohnowang(X,parameter,ink_flag,ntens,ndi,CEL)
% Funktion berechnet nachgiebigkeit im elastisch-plastischen

% Abbildungen
if ntens == 6
    [P,Paa] = set_maps(ntens,ndi);
    Pna = diag(ones(1,6));
elseif ntens == 3
    [P,Paa,~,~,Pna] = set_maps(ntens,ndi);
else
    [P,Paa,Pna] = set_maps(ntens,ndi);
end
% Parameter
M = (length(parameter)-3)/3;                                               % Anzahl TBST
r0 = parameter(end);                                                       % startradius fliessfläche
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);
% Zustandsvariablen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
    % epsp = X(ntens+1:2*ntens);
end
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% Normale an die Fließfläche
n = sqrt(3/2) *diag([1,1,2])* (s - a) / r0;
% tangentenmodul
h = tangmod(c_i,r_i,chi_i,alpha,n,Pna,Paa);
% Cel angewendet auf normale
dummy = CEL * n;
dummy2 = n' * dummy;
% Tangetiale Nachgiebigkeit
CEP = CEL - (dummy * dummy')/(h+dummy2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           tangentenmodul                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = tangmod(c_i,r_i,chi_i,alphai,n,Pna,Paa)
% berechnet plastisches Tangentenmodul
h = sum(c_i.*r_i);
M = length(r_i);
norm_ai = sqrt( sum( (Paa*alphai) .* alphai) );
norm_ai(norm_ai == 0) = 1e-40;
transn = n.';
Li = alphai./norm_ai;
var = norm_ai./r_i;
var4 = transn * Pna * Li;
var4 = 0.5 * (var4 + abs(var4));
var2 = transn * Pna * alphai;
for ii = 1 : M
    h = h - ...
       c_i(ii) * var(ii)^chi_i(ii) * var4(ii) * var2(ii);  
end

end
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 3D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          materialmodell3d(~, X, ink, ink_flag, parameter , C, D, ...
                              P, P_hat, P_tilde)
% Konkretes Materialmodell für 3D Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach Ohno Wang
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
M = (length(parameter)-3)/3;                                               % Anzahl TBST

% elastische 

r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);



% oft verwendete Konstanden
w3d2 = sqrt(3/2);
delta = 1e-40;

% -------------------------------------------------------------------------
%                        Zustandsvariablen
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
%                   Normale an die Fließfläche
% -------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r0).*P_hat*(s-a);
transn = n.';

% -------------------------------------------------------------------------
%                   Ableitung der Teilbacksztresstensoren
% -------------------------------------------------------------------------

% normen der Teilbackstresstensoren 
norm_ai = sqrt( sum( (P_hat*alpha) .* alpha) );
norm_ai(norm_ai == 0) = delta;

% % Hilfsgrößen
% Li = alpha./norm_ai;
% dummy1 = norm_ai./r_i;
% dummy1( dummy1 > 1 ) = 1;
% dummy1 = dummy1.^chi_i;
% dummy2 = transn * Li;
% dummy2 = 0.5 * ( dummy2 + abs(dummy2));
% dummy3 = P_tilde * n * r_i;
% % Ableitung Teilbackstresstensoren
% dalpha_dp = c_i .* ( dummy3 - dummy1 .* dummy2 .* alpha );

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
%                   plastischer Tangentenmodul
% -------------------------------------------------------------------------

h = transn * da_dp ;

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

% -------------------------------------------------------------------------
%                   zusammenfassen der Inkremente
% -------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dp];   

end % Ende Modell 3D


























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen Ebener Spannungszustand             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          materialmodellESZ(~, X, ink, ink_flag, parameter, C, D, P, ...
                            P_line, P_hat, A, P_check)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach OhnoWang
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
M = (length(parameter)-3)/3;                                               % anzahl Backstresstensoren

r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);

% oft verwendete Konstanden
w3d2 = sqrt(3/2);

delta = 1e-40;                                                             % numerisch 0

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
    % epsp = X(ntens+1:2*ntens);
end

% backstresstensoren
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r0).* P_hat * (s-a);
transn = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Normen der Teilbackstresstensoren
norm_ai = sqrt( sum( (P_line*alpha) .* alpha ) );
norm_ai(norm_ai == 0) = delta;                                             % Abfangen von 0
% % Hilfsgrößen
% Li = alpha./norm_ai;
% dummy1 = norm_ai./r_i;
% dummy1( dummy1 > 1 ) = 1;
% dummy1 = dummy1.^chi_i;
% dummy2 = transn * P_check * Li;
% dummy2 = 0.5 * (dummy2 + abs(dummy2));
% dummy3 = A * n * r_i;
% % Ableitungen teilbackstress
% dalpha_dp = c_i .* (dummy3 - dummy1 .* dummy2 .* alpha); 

% Hilfsvariablen 
% % Li = alpha./r_i;
% Li = alpha./norm_ai;
% h_i = c_i .* r_i;
% var2 = A * n;
% var3 = norm_ai./r_i;
% var3(var3 >= 1) = 1;
% % var3(var3 < 1 ) = 0;
% var4 = transn * P_check * Li;
% var4 = 0.5 * (var4 + abs(var4));

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

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = transn * P_check * da_dp;

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

depsp=dp.*n;
dalpha=dalpha_dp.*dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
% dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dp]; 
% dX = [dout; depsp ; dalpha([1:M*ntens]).'; dp]; 
dX = [dout; depsp ; dalpha(:); dp];

end % Ende Modell ESZ
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen Ebener Spannungszustand             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          materialmodellST(~, X, ink, ink_flag, parameter, C, D, P, ...
                            P_line, P_hat, A, P_check)
% Konkretes Materialmodell für sigma-tau Spannungszustände, gibt bei 
% vorgabe eines Lastinkrementes die Inkremente der inneren Variablen 
% zurück. Dabei wird % angenommen, das jedes übergebene Inkrement 
% elastisch-plastische Deformationen hervorruft.
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach OhnoWang
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
M = (length(parameter)-3)/3;                                               % anzahl Backstresstensoren

r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
chi_i = parameter(3+2*M:2+3*M);

% oft verwendete Konstanden
w3d2 = sqrt(3/2);

delta = 1e-40;                                                             % numerisch 0

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
    % epsp = X(ntens+1:2*ntens);
end

% backstresstensoren
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P .* sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r0).* P_hat .* (s-a);
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

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = transn * (P_check .* da_dp);

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

depsp=dp.*n;
dalpha=dalpha_dp.*dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
% dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dp]; 
% dX = [dout; depsp ; dalpha([1:M*ntens]).'; dp]; 
dX = [dout; depsp ; dalpha(:); dp];

end % Ende Modell Sigma-tau Spannungszustand

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 1D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = materialmodell1d(~, X, ink, ink_flag, parameter)
% Konkretes Materialmodell für 1D Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach OhnoWang
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

% -------------------------------------------------------------------------
%                      Materialparameter
% -------------------------------------------------------------------------

% Zuweisen der Materialparameter
M = (length(parameter)-3)/3; % Anzahl TBST
E = parameter(1);
c_i = parameter(3:M+2);
r_i = parameter(M+3:2*M+2);
chi_i = parameter(2*M+3:3*M+2);
r0 = parameter(end); % startradius fliessfläche

% -------------------------------------------------------------------------
%                      Zustandsvariablen
% -------------------------------------------------------------------------

% Auslesen der Zustandsvariable
if ink_flag == 0 % spanste
    eps = X(1);
    epsp = X(2);
    sig = E * (eps - epsp);
elseif ink_flag ==1 %dehnsteu
    sig = X(1);
end
alphai = X(3:2+M)';
alpha = sum(alphai);

% -------------------------------------------------------------------------
%                      Normale an die Fließfläche
% -------------------------------------------------------------------------

n = sqrt(3/2) * 1/r0 * (sig-alpha);

% -------------------------------------------------------------------------
%                      Ableitung Teilbackstress
% -------------------------------------------------------------------------

% normen der teilbackstresstensoren
norm_ai = sqrt(2/3) * abs(alphai);
% Hilfsvariablen
Li = sqrt(3/2) * sign(alphai);
dummy1 = norm_ai./r_i;
dummy1(dummy1 > 1) = 1;
dummy1 = dummy1.^chi_i;
dummy2 = 2/3 * n * Li;
dummy2(dummy2 < 0) = 0;
% Ableitung teilbackstresstensor
dalphaidp = c_i .* ( r_i .* n - dummy1 .* dummy2 .* alphai);
% Ableitung gesamtbackstress
dalphadp = sum(dalphaidp);

% -------------------------------------------------------------------------
%                      plastischer Tangentenmodul
% -------------------------------------------------------------------------

h = 2/3 * n * dalphadp;

% -------------------------------------------------------------------------
%                      inkrement plastische Bogenlänge
% -------------------------------------------------------------------------

if ink_flag == 0 % Spansteu
    dp = 2/3 * n * ink / h;
else % dehnsteu
    dp = 2/3 * n * E * ink / (h + 2/3 * E);
end

% -------------------------------------------------------------------------
%                      inkrement Zustandsvariablen
% -------------------------------------------------------------------------

depsp = 2/3 * dp * n;
% inversesinkrement
if ink_flag == 0 %spansteu
    deps = ink/E + depsp;
else % dehnsteu
    dsig = E * (ink - depsp);
end
% Inkrement kinematische Verfestigung
dalphai = dalphaidp * dp;

% -------------------------------------------------------------------------
%                      Zusammenfassen der Inkremente
% -------------------------------------------------------------------------

if ink_flag == 0 %spansteu
    dX = [deps;depsp;dalphai';dp];
else % dehnsteu
    dX = [dsig;depsp;dalphai';dp];
end

end % Ende Modell 1D