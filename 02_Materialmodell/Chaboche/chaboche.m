function [X_neu,dummy,CEP] = chaboche(ntens,ndi,ink, X, ink_flag, parameter)
% Materialmodell nach Chaboche mit:
% Backstresstensoren als Reihe nach Chaboche
% Einzelne Teilbackstresstensoren nach Armstrong-Frederick Kinematik
% isotrope Verfestigung asympt. gegen Grenzwert.
%
%   INPUT:
%         ntens     -> Anzahl Tensorkomponten (in Spannungen)
%         ndi       -> Anzahl Diagonalelemente
%         nsh       -> Anzahl Nebendiagonalelemente
%         ink       -> Belastungsinkrement
%         X         -> Zustandsvariablen [eps;epsp;alphai;r;p] bei spansteu
%                                        [sig;epsp;alphai;r;p] bei dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         parameter -> Materialparameter des Modells
%                      [E,nu,q,gamma,zeta_i,r_i,r0]
%
%
%    OUTPUT:
%        X_neu -> neue Zustandsvariablen nach Lastinkrement
%    
%    NOTATIONEN:
%    Spannungen                 Verzerrungen
%
%           sig_11                     eps_11
%           sig_22                     eps_22
%     sig = sig_33              eps =  eps_33
%           sig_12                    2eps_12
%           sig_13                    2eps_13
%           sig_23                    2eps_23
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: Januar 2020                                                  |
%  ----------------------------------------------------------------------

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------
% % Windows ausgabe
% fprintf('Zustandsvariablen Input:\n')
% fprintf('%.32d \n',X)
% elastizitätsparameter
E = parameter(1);
nu = parameter(2);

%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% statische Matrizen je nach Spannungszustand
[P, P_line] = set_maps(ntens,ndi);

%--------------------------------------------------------------------------
%                  Zustandsvariablen auslesen                                 
%--------------------------------------------------------------------------

% Anzahl Backstress
M = (length(parameter)-5)/2;
% Spannungen und Dehnungen
if ink_flag == 0 % spansteu
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    C = elast_steifigkeit(E,nu,ntens,ndi);
    sig = C * (eps - epsp);
    dsig_tr = ink;    % elastischer trial step
elseif ink_flag == 1 % dehnsteu
    sig = X(1:ntens);
    C = elast_steifigkeit(E,nu,ntens,ndi);
    dsig_tr = C * ink; % elastischer trial step
else
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end
% backstresstensor
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
% radius fließfläche
r = X((M+2)*ntens+1);
    
%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------

% Trial Inkrement der Spannungen 
if ntens == 1
    s = P * sig;                                                           % Span dev
    ds = P * dsig_tr;                                                      % dev des Spannungsinkrements
    s_tr = s + ds;                                                         % Trial Spannungsdeviator
    a = sum(alpha,2);                                                      % Backstress;
    beta = s_tr - a;                                                       % Relativ Spannung
    F_tr = abs(beta) - r;
elseif ntens == 2 % sigma-tau
    s = P .* sig;                                                           % Span dev
    ds = P .* dsig_tr;                                                      % dev des Spannungsinkrements
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
        eps = eps + D*ink;
        X_neu = [eps;X(ntens+1:end)];
    elseif ink_flag == 1 % dehnsteu
        sig = sig + dsig_tr;
        X_neu = [sig;X(ntens+1:end)];
    end
    
    if nargout == 2 % berechne Elastisch plastische steifigkeitsmatrix
        dummy = 1;
    elseif nargout == 3
        dummy = 1;
        CEP = elast_steifigkeit(E,nu,ntens,ndi);
    end

%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%--------------------------------------------------------------------------

else % Ermittle elastischen Anteil xel
    
    % 1. ermittle elastischen Anteil
%     xel = elastink(s,a,r,ds,ndi);
    xel = elastink2(s,a,P_line,r,ds,FTOL);
    
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
        
        % elast nachgiebigkeit
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % Abbildungen
        [P, P_hat, P_tilde] = set_maps(ntens,ndi);
        % integration
        
        % Matlab Version explizit
%         [~,X_neu] = rk87(@materialmodell3d,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D, P, P_hat, P_tilde);
                    
        % Mex Version explizit
        [~,X_neu] = rk87(@CoderVersion_3D_Chaboche_mex,[0,1], X, options,...
                         M, ink, ink_flag, parameter, C, D, P, P_hat, P_tilde);
                    
    elseif ntens == 3 && ndi == 2 % ESZ
        
        % elast nachgiebigkeit
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % Abbildungen
        [P, ~, P_hat, A, P_check] = set_maps(ntens,ndi);
        % integration
        
        % Matlab Version
%         [~,X_neu] = rk87(@materialmodellESZ,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D, P, P_hat, A, P_check);

        % Mex Version
        [~,X_neu] = rk87(@CoderVersion_ESZ_Chaboche_mex,[0,1], X, options,...
                        M, ink, ink_flag, parameter, C, D, P, P_hat, A, P_check);
    
    elseif ntens == 2 && ndi == 1 % sigma-tau
        
        % elast nachgiebigkeit
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % Abbildungen
        [P, ~, P_hat, A, P_check] = set_maps(ntens,ndi);
        % integration
        
        % Matlab Version
%         [~,X_neu] = rk87(@materialmodellST,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D, P, P_hat, A, P_check);
        
                    
        % Mex Version
        [~,X_neu] = rk87(@CoderVersion_Chaboche_ST_mex,[0,1], X, options,...
                        M, ink, ink_flag, parameter, C, D, P, P_hat, A, P_check);            
                    
                    
    elseif ntens == 1 % 1D
        
        D = 1/C;
        
        [~,X_neu] = rk87(@materialmodell1d,[0,1], X, options, ink,...
                        ink_flag, parameter);
                    
    end
    
    % nur letzten Zeitschritt ausgeben
    X_neu = X_neu(end,:)'; 
    
    % Testen der Konsistenzbedingung
    X_neu = konsistenzbedingung(X_neu,M,P,P_line,ntens,C,D,ink_flag);
    
    
    % Falls tangentiale Steifigkeit gebraucht wird
    if nargout == 2 % berechne Elastisch plastische steifigkeitsmatrix
        dummy = 1;
    elseif nargout == 3
        dummy = 1;
        CEP = tangsteifigkeit_chaboche(X_neu,parameter,ink_flag,ntens,ndi);
    end
    
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Prüfe die konsistenzbedingung                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = konsistenzbedingung(X,M,P,P_line,ntens,C,D,ink_flag)
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
r = X(end-1);

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
















% -------------------------------------------------------------------------
% DGL des Chaboche Modells
% -------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Modellgleichungen 3D Spannungszustand (explizite Integration)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          materialmodell3d(~, X, ink, ink_flag, parameter, C, D, ...
                           P, P_hat, P_tilde)
% Konkretes Materialmodell für 3D Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach Armstrong Frederick Kinematik
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%                      = 2 + num_alpha * 2 + 1
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_hat ...   -> Diverse Abbildungen
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
M = (length(parameter)-5)/2;                                               % anzahl Backstresstensoren

% isotrope Verfestigung
q = parameter(3);
gamma = parameter(4);

% kinematische Verfestigung
zeta_i = parameter(5:4+M);
r_i = parameter(5+M:end-1);
r0 = parameter(end);                                                       % startradius fliessfläche
h_i = zeta_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);

% -------------------------------------------------------------------------
%                        Zustandsvariablen
% -------------------------------------------------------------------------

% auslesen spannungs und dehnungs zustände
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
    % epsp = X(ntens+1:2*ntens);
end
% Backstresstensoren
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
% Radius Fließfläche
r = X(end-1);

% -------------------------------------------------------------------------
%                   Normale an die Fließfläche
% -------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r).*P_hat*(s-a);
nTrans = n.';

% -------------------------------------------------------------------------
%                   Ableitung der Teilbacksztresstensoren
% -------------------------------------------------------------------------

% Ableitungen der teilbackstresstensoren
dalpha_dp = w2d3 .* P_tilde * n * h_i - zeta_i .* alpha;
% Ableitung gesamtbackstresstensor
da_dp = sum(dalpha_dp,2);

% -------------------------------------------------------------------------
%                   Ableitung des Radius
% -------------------------------------------------------------------------

dr_dp = (q-gamma*(r-r0));

% -------------------------------------------------------------------------
%                   plastischer Tangentenmodul
% -------------------------------------------------------------------------

h = w3d2 * nTrans * da_dp + dr_dp;

% -------------------------------------------------------------------------
%                   Inkrement plastische Bogenlänge
% -------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (w3d2/h) * (nTrans * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = (w3d2 .* nTrans * (C * ink)) / ...
         (h + 3/2.* nTrans * ( C * n));
else % Fehler
    msg = 'keine oder flasche angabe der Laststeuerung';
    error(msg)
end

% -------------------------------------------------------------------------
%                   Inkremente Zustandsvariableb
% -------------------------------------------------------------------------

depsp=dp.*w3d2.*n;
dalpha=dalpha_dp.*dp;
dr=dr_dp * dp;

% -------------------------------------------------------------------------
%                   zusammenfassen der Inkremente
% -------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end

dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dr; dp];     

end % Ende Modell 3D

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Modellgleichungen Ebener Spannungszustand (explizite Integration)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          materialmodellESZ(~, X, ink, ink_flag, parameter, C, D, ...
                            P, P_hat, A, P_check)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach Armstrong Frederick Kinematik
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%                      = 2 + num_alpha * 2 + 1
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

ntens = 3;                                                                 % Tensorkomponenten
M = (length(parameter)-5)/2;                                               % anzahl Backstresstensoren

% isotrope Verfestigung
q = parameter(3);
gamma = parameter(4);
% kinematische Verfestigung
zeta_i = parameter(5:4+M);
r_i = parameter(5+M:end-1);
r0 = parameter(end);                                                       % startradius fliessfläche
h_i = zeta_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% spannungen und dehnungen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
end
% backstress
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
% radius fließfläche
r = X(end-1);

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r).* P_hat * (s-a);
nTrans = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Ableitung Teilbackstresstensoren
dalpha_dp = w2d3 .* A * n * h_i - zeta_i.*alpha;
% Ableitung gesamtbackstresstensor
da_dp = sum(dalpha_dp,2);

%--------------------------------------------------------------------------
%                   Ableitung des radius der Fließfläche
%--------------------------------------------------------------------------

dr_dp = (q-gamma*(r-r0));

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = w3d2 * nTrans * P_check * da_dp + dr_dp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (w3d2/h) * (nTrans * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = (w3d2 .* nTrans * (C * ink)) / ...
         (h + 3/2.* nTrans * ( C * n));
else % Fehler
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp=dp.*w3d2.*n;
dalpha=dalpha_dp.*dp;
dr=dr_dp * dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dr; dp]; 
    
end % Ende Modell ESZ








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Modellgleichungen Ebener Spannungszustand (explizite Integration)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          materialmodellST(~, X, ink, ink_flag, parameter, C, D, ...
                            P, P_hat, A, P_check)
% Konkretes Materialmodell für sigma-tau Spannungszustände, gibt bei 
% vorgabe eines Lastinkrementes die Inkremente der inneren Variablen 
% zurück. Dabei wird angenommen, das jedes übergebene Inkrement
% elastisch-plastische Deformationen hervorruft.
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach Armstrong Frederick Kinematik
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%                      = 2 + num_alpha * 2 + 1
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2021 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

ntens = 2;                                                                 % Tensorkomponenten
M = (length(parameter)-5)/2;                                               % anzahl Backstresstensoren

% isotrope Verfestigung
q = parameter(3);
gamma = parameter(4);
% kinematische Verfestigung
zeta_i = parameter(5:4+M);
r_i = parameter(5+M:end-1);
r0 = parameter(end);                                                       % startradius fliessfläche
h_i = zeta_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% spannungen und dehnungen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
end
% backstress
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
% radius fließfläche
r = X(end-1);

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P .* sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r).* P_hat .* (s-a);
nTrans = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Ableitung Teilbackstresstensoren
dalpha_dp = w2d3 .* A .* n * h_i - zeta_i.*alpha;
% Ableitung gesamtbackstresstensor
da_dp = sum(dalpha_dp,2);

%--------------------------------------------------------------------------
%                   Ableitung des radius der Fließfläche
%--------------------------------------------------------------------------

dr_dp = (q-gamma*(r-r0));

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = w3d2 * nTrans * (P_check .* da_dp) + dr_dp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (w3d2/h) * (nTrans * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = (w3d2 .* nTrans * (C * ink)) / ...
         (h + 3/2.* nTrans * ( C * n));
else % Fehler
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp=dp.*w3d2.*n;
dalpha=dalpha_dp.*dp;
dr=dr_dp * dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dr; dp]; 
    
end % Ende Modell sigma-tau









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 1D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = materialmodell1d(~, X, ink, ink_flag, parameter)
% Konkretes Materialmodell für 1D Spannungszustände, gibt bei vorgabe eines
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
%
%   OUTPUT:
%         dX        -> Inkrement der Zustandsvariablen
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%                      Materialparameter
% -------------------------------------------------------------------------

% Zuweisen der Materialparameter
M = (length(parameter)-5)/2;                                               % anzahl Backstresstensoren
% elastische 
E = parameter(1);                                                          % E-Modul
% isotrope Verfestigung
q = parameter(3);
gamma = parameter(4);
% kinematische Verfestigung
zeta_i = parameter(5:4+M);
r_i = parameter(5+M:end-1);
r0 = parameter(end);                                                       % startradius fliessfläche
h_i = zeta_i.*r_i;

% -------------------------------------------------------------------------
%                      Zustandsvariablen
% -------------------------------------------------------------------------

% Auslesen der Spannungen und Dehnungen
if ink_flag == 0 % spanste
    eps = X(1);
    epsp = X(2);
    sig = E * (eps - epsp);
elseif ink_flag ==1 %dehnsteu
    sig = X(1);
end
% Backstresstensoren
alphai = X(3:2+M)';                                                        % einzelne Anteile
alpha = sum(alphai);                                                       % Gesamter Backstresstensor
% radius Fließfläche
r = X(3+M);

% -------------------------------------------------------------------------
%                      Ableitung Radius Fließfläche
% -------------------------------------------------------------------------

drdp = (q-gamma*(r-r0));

% -------------------------------------------------------------------------
%                      Ableitung Teilbackstress
% -------------------------------------------------------------------------

% Ableitung Teilbackstress
dalphaidp = (h_i .* (sign(sig-alpha) - 1./r_i .* alphai)).';
% Ableitung Gesamtbackstress
dalphadp = sum(dalphaidp);

% -------------------------------------------------------------------------
%                      plastischer Tangentenmodul
% -------------------------------------------------------------------------

h = sign(sig-alpha) * dalphadp + drdp;

% -------------------------------------------------------------------------
%                      inkrement plastische Bogenlänge
% -------------------------------------------------------------------------

if ink_flag == 0 % Spansteu
    dp = sign(sig-alpha) * ink / h;
elseif ink_flag == 1 % dehnsteu
    dp = sign(sig - alpha) * E * ink / (h + E);
end
% -------------------------------------------------------------------------
%                      inkrement Zustandsvariablen
% -------------------------------------------------------------------------

% plastische dehnung
depsp = dp * sign(sig - alpha);
% spannungen oder dehnungen
if ink_flag == 0 %spansteu
    deps = ink/E + depsp;
elseif ink_flag == 1 % dehnsteu
    dsig = E * (ink - depsp);
end
% Inkrement isotrope verfestigung
dr = drdp * dp;
% Inkrement kinematische Verfestigung
dalphai = dalphaidp * dp;

% -------------------------------------------------------------------------
%                      Zusammenfassen der Inkremente
% -------------------------------------------------------------------------

if ink_flag == 0 %spansteu
    dX = [deps;depsp;dalphai;dr;dp];
else % dehnsteu
    dX = [dsig;depsp;dalphai;dr;dp];
end
end % Ende Modell 1D


























% -------------------------------------------------------------------------
% Hilfsfunktionen
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Tangentialer elastisch-plastischer Steifigkeitstensor
function CEP = tangsteifigkeit_chaboche(X_neu,parameter,ink_flag,ntens,ndi)
        % Matrizen je nach Spannungszustand
        P = set_maps(ntens,ndi);
        % Parameter 
        M = (length(parameter)-5)/2;
        E = parameter(1);
        nu = parameter(2);
        q = parameter(3);
        gamma = parameter(4);
        % kinematische Verfestigung
        zeta_i = parameter(5:4+M);
        r_i = parameter(5+M:end-1);
        r0 = parameter(end); % startradius fliessfläche
        h_i = zeta_i.*r_i;
        CE = elast_steifigkeit(E,nu,ntens,ndi);
        % Neue Zustände
        if ink_flag == 0 % spansteu
            eps = X_neu(1:ntens);
            epsp = X_neu(ntens+1:2*ntens);
            sig = CE * (eps - epsp);
        elseif ink_flag == 1 % dehnsteu
            sig = X_neu(1:ntens);
        end
        alpha = zeros(ntens,M);
        for i = 1:M
            alpha(:,i) = X_neu((i+1)*ntens+1:(i+2)*ntens);
        end
        r = X_neu((M+2)*ntens+1);
        % berechne normale
        s = P * sig;
        a = sum(alpha,2);
        beta = s - a;
        n = (sqrt(3/2)/r).*diag([ones(1,ndi),2*ones(1,ntens-ndi)])*beta;
        % plastischer Modul 
        h = tangmod_chaboche(n,alpha,r,r0,h_i,zeta_i,q,gamma,ntens,ndi);
        % tangentiale Steifigkeit
        dummy = CE * n;
        CP = 3/2 .*(dummy * dummy')/(3/2 * n' * dummy + h) ;
        CEP = CE - CP;
end

% -------------------------------------------------------------------------
% Tangenten Modul Chaboche , berechnet den Tangentenmodul fürs Chaboche
% Modell
function h = tangmod_chaboche(n,alphai,r,r0,h_i,zeta_i,q,gamma,ntens,ndi)
% Anzahl Backstresstensoren
M = length(h_i);
% Spannungsdev je nach Spannungszustand
if ntens == 1
    P_line = 1;
elseif ntens == 3 && ndi == 2
    P_line = [2,1,0;...
             1,2,0;...
             0,0,1];
elseif ntens == 6
    P_line = diag([1,1,1,1,1,1]);
end   
h = sum(h_i) + q - gamma * (r-r0);
for i = 1:M
    h = h - sqrt(3/2) * zeta_i(i) * n' * P_line * alphai(:,i);
end
end

