function [Xneu,CEP,DEP] = doring(ntens, ndi, ink, X, ink_flag, para)
% Materialmodell von Döring
%
% Quelle:
%       Zum Deformation und Schädigungsverhalten metallischer Werkstoffe
%       unter mehrachsigen nichtproportionalen zyklischen Beanspruchungen
%       - Dissertation Döring 2006
%
%
%   INPUT:
%         ntens     -> Anzahl Tensorkomponten (in Spannungen)
%         ndi       -> Anzahl Diagonalelemente, mehr kkommentare
%         ink       -> Belastungsinkrement
%         X         -> Zustandsvariablen 
%         X = [eps;epsp;alphai;beta;p;rm;A;ri;CT] bein spansteu
%         X = [sig;epsp;alphai;beta;p;rm;A;ri;CT] bein dehnsteu
%         X(1:ntens)   = eps                                                   -> Gesamtdehnung
%         X(1:ntens)   = sig                                                   -> Spannungen
%         X(ntens+1:2*ntens)  = epsp                                           -> plastische Dehnungen
%         X(2*ntens+1:(2+M)*ntens) = alpha_i                                   -> Teilbackstresstensoren
%         X((2+M)*ntens+1:(3+M)*ntens) = beta                                  -> (Backstraintensor) Mittelpunkt Gedächtnisfläche im Dehnungsraum
%         X((3+M)*ntens+1) =  p                                                -> plastische Bogenlänge, hier dp = sqrt(dep:dep)
%         X((3+M)*ntens+2) = rm                                                -> Dehnungsgedächtnisfläche
%         X((3+M)*ntens+3) = A                                                 -> Nichtproportionaalitätskennwert
%         X((3+M)*ntens+4) = sigF                                              -> Fließspannung, r0 * sqrt(3/2) = sig_F
%         X((3+M)*ntens+5 : (3+M)*(ntens+1) + 1) = ri                          -> Begrenzungsradien für Teilbackstresstensoren
%         X((3+M)*(ntens+1) + 2 : (3+M+ntens/2)*(ntens+1) + 1) = CT            -> Nichtproportionalitätstensor nach Tanaka
%                                                                             
%
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         para      -> Materialparameter des Modells
%                      siehe Quelle (Seite 105 ff)
%         ... Elastizität
%         para(1) = E           -> E-Modul
%         para(2) = nu          -> Querdehnzahl
%         ... zyklisch stabilisiert
%         para(3:2+M) = ci         (i=1...M) 
%         para(3+M:3+2*M) = rinfi  (i=0...M)
%         ... transiente zyklische Ver und Entfestigung
%         para(4+2*M:4+3*M) = a1i  (i=0...M)
%         para(5+3*M:5+4*M) = a2i  (i=0...M)
%         para(6+4*M:6+5*M) = a3i  (i=0...M)
%         para(7+5*M) = b1  
%         para(8+5*M) = b2  
%         para(9+5*M) = b3 
%         para(10+5*M) = r0init
%         ... nichtproportionale Zusatzverfestigung
%         para(11+5*M) = qn0
%         para(12+5*M) = an
%         para(13+5*M) = bn
%         para(14+5*M) = ct
%         para(15+5*M) = ca
%         ... Non Masing Verhalten
%         para(16+5*M:16+6*M) = api (i=0...M)
%         para(17+6*M) = bp
%         ... verhalten der Dehnungsgedächtnissfläche
%         para(18+6*M) = eta
%         para(19+6*M) = cme
%         para(20+6*M) = omega
%         ... proportionales Ratchetting
%         para(21+6*M:20+7*M) = Qi  (i=1...M)
%         para(21+7*M) = achi
%         para(22+7*M) = bchi
%         ... nichtproportionales Ratchetting
%         para(23+7*M:22+8*M) = cchii      (i=1...M)
%         
%
%    OUTPUT:
%        X_neu -> neue zustandsvariablen nach Lastinkrement
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
%
%    Tanaka Tensor
%     -> Symmetrischer Tensor 4. Stufe (Dehnungscharakter)
%     -> 21 unabhängige Komponenten (bei ntens = 6)
%     -> Speicherrichtung entlang der Diagonalen
%            Tensorkomp. in Voigt Notation                      Für Zustandsvariablen
%        ( 1111   1122   1133  2*1112  2*1113  2*1123 )       ( 1  7  12  16  19  21)
%        (        2222   2233  2*2212  2*2213  2*2223 )       (    2   8  13  17  20)
%   CT = (               3333  2*3312  2*3312  2*3323 )  =    (        3   9  14  18)
%        (                     4*1212  4*1213  4*1223 )       (            4  10  15)
%        (                             4*1313  4*1323 )       (                5  11)
%        (                                     4*2323 )       (                    6)
%
%
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: September 2020                                               |
%  ----------------------------------------------------------------------



%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------
M = (length(para)-22)/8;                                                   % Anzahl Backstress 
E = para(1);                                                               % E-Modul
nu = para(2);                                                              % Querdehnzahl
C = elast_steifigkeit(E,nu,ntens,ndi);                                     % Elastische Steifigkeitsmatrix

%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% statische Matrizen
[P, PLINE] = set_maps(ntens,ndi);


%--------------------------------------------------------------------------
%            Identifikation der Zustandsvariablen                         
%--------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteuerung
    eps = X(1:ntens);                                                      % Dehnung
    epsp = X(ntens+1:2*ntens);                                             % plast. Dehnung
    sig = C *(eps-epsp);                                                   % Spannung
    dsig_tr = ink;                                                         % (elast.trial)Spannungsinkrement
else             % Dehnungssteurung
    sig = X(1:ntens);                                                      % Spannugen
    dsig_tr = C*ink;                                                       % (elast.trial)Spannungsinkrement
end

alpha = reshape(X(2*ntens+1:(2+M)*ntens),ntens,M);                         % Backstresstensoren

%--------------------------------------------------------------------------
%               Radius der Fließfläche berechnen                          
%--------------------------------------------------------------------------
r0 = X((3+M)*ntens+4);

%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------
s = P*sig;                                                                 % Spannungsdeviator
a = sum(alpha,2);                                                          % Backstresstensor
ds = P * dsig_tr;                                                          % (trial)Dev.Span.Ink
s_tr = s + ds;                                                             % (trial)Deviatorspannungen
eff = s_tr - a;                                                            % (trial)Effektivspannung
F_tr = eff'*PLINE*eff - 2/3 * r0^2;                                               % (trial) Fließfunktion
FTOL = 1e-7;                                                               % Toleranz für Verletzung der Konsitenzbedinung

%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------
if F_tr < FTOL % Trial Step annehmen
    % ... elastisches update der zustandsvariablen
    D = elast_nachgiebigkeit(E,nu,ntens,ndi);
    if ink_flag == 0 % Spannsteuerung
        eps = eps + D * ink;
        Xneu = [eps;X(ntens+1:end)];
    else             % Dehnungssteuerung
        sig = sig + dsig_tr;
        Xneu = [sig;X(ntens+1:end)];
    end
    
    % ... Ausgabe Steifigkeiten
    if nargout == 2
        CEP = C;
    elseif nargout == 3
        CEP = C;
        DEP = D;
    end
%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%-------------------------------------------------------------------------- 
else % Trial step ablehnen
    % ... elastischen Anteil ermitteln
%     xel = elastink(s,a,abs(r0),ds,ndi);
    xel = elastink2(s,a,PLINE,abs(r0),ds,FTOL);
    % ... elastischen Anteil aufbringen
    if ink_flag == 0 % Spannsteuerung
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        eps = eps + D * (xel*ink);
        X = [eps;X(ntens+1:end)];
    else             % Dehnungssteuerung
        sig = sig + xel * dsig_tr;
        X = [sig;X(ntens+1:end)];
    end
    % ... plastischer Anteil des Inkrements
    ink = (1-xel)*ink;
    % ... integration je nach Spannungszustand
    options = [];                                                          % optionen für Runge Kutta Verfahren
    
    % ... Ändere Fließspannung auf dev. Fließspannung
    X((3+M)*ntens+4) = X((3+M)*ntens+4) * sqrt(2/3);
    if ntens == 6 % 3D
        % ... elastische Nachgebigkeit
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % ... Abblidungen
        [P, PLINE, PHAT] = set_maps(ntens,ndi);
        PCHECK = diag([ 1 1 1 1 1 1]);
        AHAT = PLINE;
        % ... Integration
        [~,Xneu] = rk87(@doring3D,[0,1], X, options, ink,...
                   ink_flag, para, C, D, P, PLINE, PHAT);
    elseif ntens == 3 && ndi == 2 % ESZ
        % ... elastische Nachgebigkeit
        D = elast_nachgiebigkeit(E,nu,ntens,ndi);
        % ... Abblidungen
        [P, PLINE, PHAT, AHAT, PCHECK] = set_maps(ntens,ndi);
        % ... Integration
        [~,Xneu] = rk87(@doringESZ,[0,1], X, options, ink,...
                   ink_flag, para, C, D, P, PLINE, PHAT, AHAT, PCHECK);
    else
        msg = 'Spannungszustand nicht implementiert';
        error(msg)
    end % Ende Verzweigung Spannungszustände
    
    % ... nur letzten Schritt ausgeben
    Xneu = Xneu(end,:)'; 
    
    % ... Ausgabe Steifigkeiten
    if nargout == 2
        CEP = TangSteifEP(ntens,M,ink_flag,Xneu,para,...
                           C,P,PHAT,PLINE,PCHECK,AHAT);
    elseif nargout == 3
        CEP = TangSteifEP(ntens,M,ink_flag,Xneu,para,...
                           C,P,PHAT,PLINE,PCHECK,AHAT);
        DEP = TangNachEP(ntens,M,ink_flag,Xneu,para,...
                           D,P,PHAT,PLINE,PCHECK,AHAT);
    end
    
    % ... Ändere Deviatorische FS wieder auf nichtdev. FS
    Xneu((3+M)*ntens+4) = Xneu((3+M)*ntens+4) * sqrt(3/2);
    
    
    
end % Ende Verzweigung elastischer Trial Step
end % Ende der Hauptfunktion




















































% ----------------------------------------------------------------------- %
%                Modellgleichungen 3D Spannungszustand                    %
% ----------------------------------------------------------------------- %
function dX = doring3D(~,X,ink,ink_flag,para,...
                       C,D,P,PLINE,PHAT)
% Funktion berechnet Inkremente der Zustandsvariablen fürs Döringmodell im
% 3D Spannungszustand
% 
% INPUT:
%     t          - Zeit(hier nicht gebraucht nur für ode solver)
%     X          - Zustandsvariablen
%    ink         - Lastinkrement
%   ink_flag     - Laststeuerung
%    para        - Modellparameter (siehe oben)
%     C          - elastische Steifigkeit
%     D          - elastische Nachgebigkeit
%     P          - Abbildung von Spannungsraum in dev.Spannungsraum
%   PLINE        - Skalarproduk von Deviatorspannungen
%   PHAT         - Abbildung Dehnungen zu Spannungen
%   .....
% 
% OUTPUT:
% dX      - Inkremente aller Zustandsvariablen X in der oben festgelegten 
%           Reihenfolge
%__________________________________________________________________________

% -------------------------------------------------------------------------
% Identifikation der Modellparameter
% -------------------------------------------------------------------------
ntens = 6;                                                                 % Anzahl Tensorkomponenten
M = (length(para)-22)/8;                                                   % Anzahl Backstress 
% ... Elast.
E = para(1);
nu = para(2);
% ... zyk.stab.
ci = para(3:2+M);
rinfi = para(3+M:3+2*M);
% ... transiente zyk. Verf.
a1i = para(4+2*M:4+3*M);
a2i = para(5+3*M:5+4*M);
a3i = para(6+4*M:6+5*M);
b1 = para(7+5*M);
b2 = para(8+5*M);
b3 = para(9+5*M);
r0init = para(10+5*M);
% ... nichtproportionale Zusatzverfestigung
qn0 = para(11+5*M);
an = para(12+5*M);
bn = para(13+5*M);
ct = para(14+5*M);
ca = para(15+5*M);
% ... non Masing Verhalten
api = para(16+5*M:16+6*M);
bp = para(17+6*M);
% ... Dehnungsgedächtnissfläche
eta = para(18+6*M);
cme = para(19+6*M);
omega = para(20+6*M);
% ... proportionales Ratchetting
Qi = para(21+6*M:20+7*M);
achi = para(21+7*M);
bchi = para(22+7*M);
% ... nichtproportionales Ratchetting
cchii = para(23+7*M:22+8*M);
% ... default Korrekturterme
[aki,bk] = defaultak(M);
b = 200; % default geschwindigkeit für annäherung der radien an die Target Werte


% -------------------------------------------------------------------------
% Identifikation der Zustandsvariablen
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteurung
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps-epsp);
else             % Dehnungssteurung
    sig = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
end
% ... Backstress
alpha = reshape(X(2*ntens+1:(2+M)*ntens),ntens,M);
% ... Backstrain
beta = X((2+M)*ntens+1:(3+M)*ntens);
% ... akk. plast. Dehn
p = X((3+M)*ntens+1);
% ... Gedäch.Fläche Dehn.raum
rm = X((3+M)*ntens+2);
% ... NP Para
A = X((3+M)*ntens+3);
% ... Radius FF in dev. Spannungsraum
r0 = X((3+M)*ntens+4);
% ... Grenzradien Backstresstensoren in dev. Spanungsraum
ri = X((3+M)*ntens+5 : (3+M)*(ntens+1) + 1);
% ... Tanaka Tensors ( als Vektor gespeichert )
CT = X((3+M)*(ntens+1) + 2 : (3+M+ntens/2)*(ntens+1) + 1);

% -------------------------------------------------------------------------
% Backstresstensoren
% -------------------------------------------------------------------------
delta = 1e-40; % numerische Null
% ... abfangen exakte null bei radien
ri( ri==0 ) = delta;
% if rm == 0
%     rm = delta;
% end
% ... Gesamtbackstresstensor
a = sum(alpha,2);           
% norma = sqrt( sum( (PLINE*a).*a ));
% if norma == 0
%     norma = delta;
% end
% ... Normen Teilbackstresstensoren
normai = sqrt( sum( (PLINE*alpha).*alpha ) );
normai(normai == 0) = delta;
% ... Normierte Teilbackstresstensoren
Li = alpha./normai;

% -------------------------------------------------------------------------
% Fließflächennormale
% -------------------------------------------------------------------------
% ... Spannungsdeviator
s = P*sig;
% ... effektivspannung
eff = s - a;
normeff = sqrt( sum( (PLINE*eff).*eff ));
% ... Fließflächennormale (als Dehnung mit dopplten Nebendiagonalelementen)
n = PLINE* eff./normeff;
nTrans = n'; % Transposition weil die Später dauernt gebraucht wird

% -------------------------------------------------------------------------
% Berechne Exponenten chi
% -------------------------------------------------------------------------
% ... chi0
chii0 = Qi .* ( 1 + ( achi/(1+bchi*rm)^2) );
% ... chi
chii = chii0 + (chii0+0.1).*(cchii-1).*(1-abs(nTrans*Li));
chii( chii <= 0 ) = 0;

% -------------------------------------------------------------------------
% Parameterfunktion/ Target Values der Radien
% -------------------------------------------------------------------------
% ... nichtproportional
qn = qn0 * (1+an/(1+bn*rm)^2);
% ... proportional
qpi = 1+api./(1+bp*rm)^2;
% ... Endwert radius Fließfläche
rinfi_0 = rinfi(1) * ( A*qn + (1-A)*qpi(1) );
% ... Endwert radius gedächtnissfläche
rinfi_i = rinfi(2:end) .* qpi(2:end) .* ( 1 + aki./(1+bk*chii) );          % So stehts im Code
% rinfi_i = rinfi(2:end) .* qpi(2:end) .* ( 1 + aki./(1+bk*chii).^2 );          % So stehts in der Diss.
% ... Targetvalues
dummy =  1 + a1i./(1+b1*p)^2 + a2i./(1+b2*p)^2 + a3i./(1+b3*p)^2;
rT0 = rinfi_0 * dummy(1);
rTi = rinfi_i .* dummy(2:end);

% -------------------------------------------------------------------------
% Ableitungen ri i=0...M
% -------------------------------------------------------------------------
% ... Grenzeradius FF
dr0dp = rT0 - r0; 
if abs(dr0dp) <= 1
    dr0dp = b*dr0dp;
else
    dr0dp = b*abs(dr0dp)*dr0dp;
end
% ... Grenzradien Backstresstensoren
dridp = rTi' - ri;
idx = abs(dridp) <= 1;
dridp(idx) = b*dridp(idx); 
dridp(~idx)= b*abs(dridp(~idx)).*dridp(~idx);

% -------------------------------------------------------------------------
% Ableitungen Backstresstensoren
% -------------------------------------------------------------------------
% ... Wichtungsfunktionen
Wi = normai./abs(ri');
Wi(Wi > 1) = 1;
Wi = Wi.^chii;
% ... Schleife über Teilbackstresstensoren
dalphadp = zeros(ntens,M);
dummy = PHAT * n;
for i = 1:M 
    dalphadp(:,i) = ci(i) * ( ri(i)*dummy - Wi(i)*alpha(:,i) ) + ...
        alpha(:,i)./ri(i)*dridp(i);
end

% dummy = PHAT * n;
% dalphadp = ci.* ( ri.*dummy - Wi.*alpha) + alpha./ri.*dridp;

% ... gesamtbackstresstensor
dadp = sum(dalphadp,2);

% -------------------------------------------------------------------------
% plastischer Tangentenmodul
% -------------------------------------------------------------------------
h = nTrans * dadp + dr0dp;
% ... Begrenze Entfestigung
if h < -0.6 *E
    h = -0.6*E;
    dr0dp = h - nTrans * dadp;
end

% -------------------------------------------------------------------------
% Inkrement plastische Bogenlänge/akk. plastische Dehnung
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteuerung
    dp = (nTrans*ink)/h;
else             % Dehnungssteuerung
    dp = (nTrans*(C*ink))/(h + nTrans*(C*n));
end

% -------------------------------------------------------------------------
% Fließregel (inkrement plastische Dehnung)
% -------------------------------------------------------------------------
depsp = dp*n;

% -------------------------------------------------------------------------
% Spannungs/Dehnungsinkrement
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteuerung
    % ... Dehnungsinkrement
    dout = D * ink + depsp;
else             % Dehnungssteuerung
    % ... Spannungsinkrement
    dout = C * (ink - depsp);
end

% -------------------------------------------------------------------------
% Inkremente Backstresstensoren und Radien
% -------------------------------------------------------------------------
dalpha = dalphadp * dp;
dr0 = dr0dp * dp;
dri = dridp * dp;

% -------------------------------------------------------------------------
% Inkremente der Dehnungsgedächtnissfläche
% -------------------------------------------------------------------------
% ... differenz plasrische Dehnung und Mittelpunkt
dstar = (epsp - beta);
normstar2 = dstar'*PHAT*dstar;
% ... Gedächtnissfläche
Fstar = normstar2 - rm^2;
% ... Hilfsfunktionen
HF = myheaviside(Fstar);
if normstar2 == 0
    normstar2 = delta;
end
if HF == 0
    % ... ||n*|| < 1
    nstar = dstar./sqrt(normstar2);
    drm = - (cme * rm)^omega; 
else
    % ... ||n*|| = 1
    if normstar2 > delta
        nstar = dstar./sqrt(normstar2);
    else
        nstar = n;
    end
    drm = eta;
end
% ... <n*:n>
nn = nTrans * PHAT * nstar;
nn = mymacaulay(nn);

% ... inkrement gedächtnissfläche
drm = drm * nn * dp;
% ... inkrement Backstrain
dbeta = (1-eta) * HF * nn * nstar * dp;

% -------------------------------------------------------------------------
% Änderung des Tanaka Tensors und des NP Parameters
% -------------------------------------------------------------------------
[dA,dCT] = inkTanaka(ntens,A,CT,n,ca,ct,dp);

% -------------------------------------------------------------------------
% Zusammenfassen der Outputgrößen
% -------------------------------------------------------------------------
%         X = [eps;epsp;alphai;beta;p;rm;A;ri;CT] bein spansteu
dX = [dout;depsp;reshape(dalpha,ntens*M,1);dbeta;dp;drm;dA;dr0;dri;dCT];
end % Ende 3D Spannungszustand
































function dX = doringESZ(~,X,ink,ink_flag,para,...
                        C, D, P, PLINE, PHAT, AHAT, PCHECK)
% Funktion berechnet Inkremente der Zustandsvariablen fürs Döringmodell im
% ebenen Spannungszustand
% 
% INPUT:
%     t          - Zeit(hier nicht gebraucht nur für ode solver)
%     X          - Zustandsvariablen
%    ink         - Lastinkrement
%   ink_flag     - Laststeuerung
%    para        - Modellparameter (siehe oben)
%    C,D         - Elastischer Steifigkeits- & Nachgebigkeitstensor
%     P          - Abbildung Spannung auf Spannungsdeviator
%   ...
% 
% OUTPUT:
% dX      - Inkremente aller Zustandsvariablen X in der oben festgelegten 
%           Reihenfolge
%__________________________________________________________________________
% -------------------------------------------------------------------------
% Identifikation der Modellparameter
% -------------------------------------------------------------------------
ntens = 3;                                                                 % Anzahl Tensorkomponenten
M = (length(para)-22)/8;                                                   % Anzahl Backstress 
% ... Elast.
E = para(1);
nu = para(2);
% ... zyk.stab.
ci = para(3:2+M);
rinfi = para(3+M:3+2*M);
% ... transiente zyk. Verf.
a1i = para(4+2*M:4+3*M);
a2i = para(5+3*M:5+4*M);
a3i = para(6+4*M:6+5*M);
b1 = para(7+5*M);
b2 = para(8+5*M);
b3 = para(9+5*M);
r0init = para(10+5*M);
% ... nichtproportionale Zusatzverfestigung
qn0 = para(11+5*M);
an = para(12+5*M);
bn = para(13+5*M);
ct = para(14+5*M);
ca = para(15+5*M);
% ... non Masing Verhalten
api = para(16+5*M:16+6*M);
bp = para(17+6*M);
% ... Dehnungsgedächtnissfläche
eta = para(18+6*M);
cme = para(19+6*M);
omega = para(20+6*M);
% ... proportionales Ratchetting
Qi = para(21+6*M:20+7*M);
achi = para(21+7*M);
bchi = para(22+7*M);
% ... nichtproportionales Ratchetting
cchii = para(23+7*M:22+8*M);
% ... default Korrekturterme
[aki,bk] = defaultak(M);
b = 200; % default geschwindigkeit für annäherung der radien an die Target Werte


% -------------------------------------------------------------------------
% Identifikation der Zustandsvariablen
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteurung
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps-epsp);
else             % Dehnungssteurung
    sig = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
end
% ... Backstress
alpha = reshape(X(2*ntens+1:(2+M)*ntens),ntens,M);
% ... Backstrain
beta = X((2+M)*ntens+1:(3+M)*ntens);
% ... akk. plast. Dehn
p = X((3+M)*ntens+1);
% ... Gedäch.Fläche Dehn.raum
rm = X((3+M)*ntens+2);
% ... NP Para
A = X((3+M)*ntens+3);
% ... Radius FF in dev. Spannungsraum
r0 = X((3+M)*ntens+4);
% ... Grenzradien Backstresstensoren in dev. Spanungsraum
ri = X((3+M)*ntens+5 : (3+M)*(ntens+1) + 1);
% ... Tanaka Tensors ( als Vektor gespeichert )
CT = X((3+M)*(ntens+1) + 2 : (3+M+ntens/2)*(ntens+1) + 1);

% -------------------------------------------------------------------------
% Backstresstensoren
% -------------------------------------------------------------------------
delta = 1e-40; % numerische Null
% ... abfangen exakte null bei radien
ri( ri==0 ) = delta;
% ... Gesamtbackstresstensor
a = sum(alpha,2);           
% norma = sqrt( sum( (PLINE*a).*a ));
% ... Normen Teilbackstresstensoren
normai = sqrt( sum( (PLINE*alpha).*alpha ) );
normai(normai == 0) = delta;
% ... Normierte Teilbackstresstensoren
Li = alpha./normai;

% -------------------------------------------------------------------------
% Fließflächennormale
% -------------------------------------------------------------------------
% ... Spannungsdeviator
s = P*sig;
% ... effektivspannung
eff = s - a;
normeff = sqrt( sum( (PLINE*eff).*eff ));
% ... Fließflächennormale (als Dehnung mit dopplten Nebendiagonalelementen)
n = PHAT* eff./normeff;
nTrans = n'; % Transposition weil die Später dauernt gebraucht wird


% -------------------------------------------------------------------------
% Berechne Exponenten chi
% -------------------------------------------------------------------------
% ... chi0
chii0 = Qi .* ( 1 + ( achi/(1+bchi*rm)^2) );
% ... chi
chii = chii0 + (chii0+0.1).*(cchii-1).*(1-abs(nTrans*PCHECK*Li));
chii( chii <= 0 ) = 0;

% -------------------------------------------------------------------------
% Parameterfunktion/ Target Values der Radien
% -------------------------------------------------------------------------
% ... nichtproportional
qn = qn0 * (1+an/(1+bn*rm)^2);
% ... proportional
qpi = 1+api./(1+bp*rm)^2;
% ... Endwert radius Fließfläche
rinfi_0 = rinfi(1) * ( A*qn + (1-A)*qpi(1) );
% ... Endwert radius gedächtnissfläche
rinfi_i = rinfi(2:end) .* qpi(2:end) .* ( 1 + aki./(1+bk*chii) );          % So stehts im Code
% rinfi_i = rinfi(2:end) .* qpi(2:end) .* ( 1 + aki./(1+bk*chii).^2 );          % So stehts in der Diss.
% ... Targetvalues
dummy =  1 + a1i./(1+b1*p)^2 + a2i./(1+b2*p)^2 + a3i./(1+b3*p)^2;
rT0 = rinfi_0 * dummy(1);
rTi = rinfi_i .* dummy(2:end);

% -------------------------------------------------------------------------
% Ableitungen ri i=0...M
% -------------------------------------------------------------------------
% ... Grenzeradius FF
dr0dp = rT0 - r0; 
if abs(dr0dp) <= 1
    dr0dp = b*dr0dp;
else
    dr0dp = b*abs(dr0dp)*dr0dp;
end
% ... Grenzradien Backstresstensoren
dridp = rTi' - ri;
idx = abs(dridp) <= 1;
dridp(idx) = b*dridp(idx); 
dridp(~idx)= b*abs(dridp(~idx)).*dridp(~idx);

% -------------------------------------------------------------------------
% Ableitungen Backstresstensoren
% -------------------------------------------------------------------------
% ... Wichtungsfunktionen
Wi = normai./abs(ri');
Wi(Wi > 1) = 1;
Wi = Wi.^chii;
% ... Schleife über Teilbackstresstensoren
dalphadp = zeros(ntens,M);
dummy = AHAT * n;
for i = 1:M 
    dalphadp(:,i) = ci(i) * ( ri(i)*dummy - Wi(i)*alpha(:,i) ) + ...
        alpha(:,i)./ri(i)*dridp(i);
end

% dummy = PHAT * n;
% dalphadp = ci.* ( ri.*dummy - Wi.*alpha) + alpha./ri.*dridp;

% ... gesamtbackstresstensor
dadp = sum(dalphadp,2);

% -------------------------------------------------------------------------
% plastischer Tangentenmodul
% -------------------------------------------------------------------------
h = nTrans * PCHECK * dadp + dr0dp;
% ... Begrenze Entfestigung
if h < -0.6 *E
    h = -0.6*E;
    dr0dp = h - nTrans * PCHECK * dadp;
end

% -------------------------------------------------------------------------
% Inkrement plastische Bogenlänge/akk. plastische Dehnung
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteuerung
    dp = (nTrans*ink)/h;
else             % Dehnungssteuerung
    dp = (nTrans*(C*ink))/(h + nTrans*(C*n));
end

% -------------------------------------------------------------------------
% Fließregel (inkrement plastische Dehnung)
% -------------------------------------------------------------------------
depsp = dp*n;

% -------------------------------------------------------------------------
% Spannungs/Dehnungsinkrement
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteuerung
    % ... Dehnungsinkrement
    dout = D * ink + depsp;
else             % Dehnungssteuerung
    % ... Spannungsinkrement
    dout = C * (ink - depsp);
end

% -------------------------------------------------------------------------
% Inkremente Backstresstensoren und Radien
% -------------------------------------------------------------------------
dalpha = dalphadp * dp;
dr0 = dr0dp * dp;
dri = dridp * dp;

% -------------------------------------------------------------------------
% Inkremente der Dehnungsgedächtnissfläche
% -------------------------------------------------------------------------
% ... differenz plasrische Dehnung und Mittelpunkt
dstar = (epsp - beta);
normstar2 = dstar'*[2 1 0;1 2 0;0 0 0.5]*dstar;
% ... Gedächtnissfläche
Fstar = normstar2 - rm^2;
% ... Hilfsfunktionen
HF = myheaviside(Fstar);
if normstar2 == 0
    normstar2 = delta;
end
if HF == 0
    % ... ||n*|| < 1
    nstar = dstar./sqrt(normstar2);
    drm = - (cme * rm)^omega; 
else
    % ... ||n*|| = 1
    if normstar2 > delta
        nstar = dstar./sqrt(normstar2);
    else
        nstar = n;
    end
    drm = eta;
end
% ... <n*:n>
nn = nTrans * [2 1 0;1 2 0;0 0 0.5] * nstar;
nn = mymacaulay(nn);

% ... inkrement gedächtnissfläche
drm = drm * nn * dp;
% ... inkrement Backstrain
dbeta = (1-eta) * HF * nn * nstar * dp;

% -------------------------------------------------------------------------
% Änderung des Tanaka Tensors und des NP Parameters
% -------------------------------------------------------------------------
[dA,dCT] = inkTanaka(ntens,A,CT,n,ca,ct,dp);

% -------------------------------------------------------------------------
% Zusammenfassen der Outputgrößen
% -------------------------------------------------------------------------
%         X = [eps;epsp;alphai;beta;p;rm;A;ri;CT] bein spansteu
dX = [dout;depsp;reshape(dalpha,ntens*M,1);dbeta;dp;drm;dA;dr0;dri;dCT];

end % Ende Funktion Ebener Spannungszustand































% ----------------------------------------------------------------------- %
%                plastischer Tangentenmodul                               %
% ----------------------------------------------------------------------- %
function h = tangentenmodul(ntens,M,para,Ai,p,r0,ri,rm,A,n,...
                            PLINE,PCHECK,AHAT)
% Funktion berechnet plastischen Tangentenmodul aus zustandsvariablen
% INPUT:
%    ntens   - Anzahl Tentorkomponenten
%     M      - Anzahl Teilbackstresstensoren 
%  para      - Materialparameter
%
%    Ai      - TeilBackstress
%    p       - plastische Bogenlänge
%   r0       - Dev. Radius FF
%   ri       - Radien Teilbackstress
%   rm       - Dehnungsmemory
%   A        - NP Parameter
%
%   n        - Normalentensor
%
%  PLINE     - Skalarprodukt Spannungsdeviatoren
%  PCHECK    - Skalarprodukt dev Spannung, dev Dehnung (=EINS für 3D)
%  AHAT      - Dehnunen auf Spannungen (elemenieren doppelter ND
%                                       Komponenten)
%
% OUTPUT:
%   h        - plastischer Tangentenmodul
%__________________________________________________________________________
% -------------------------------------------------------------------------
% Identifikation der Modellparameter
% -------------------------------------------------------------------------
% ... Elast.
E = para(1);
% nu = para(2);
% ... zyk.stab.
ci = para(3:2+M);
rinfi = para(3+M:3+2*M);
% ... transiente zyk. Verf.
a1i = para(4+2*M:4+3*M);
a2i = para(5+3*M:5+4*M);
a3i = para(6+4*M:6+5*M);
b1 = para(7+5*M);
b2 = para(8+5*M);
b3 = para(9+5*M);
% r0init = para(10+5*M);
% ... nichtproportionale Zusatzverfestigung
qn0 = para(11+5*M);
an = para(12+5*M);
bn = para(13+5*M);
% ct = para(14+5*M);
% ca = para(15+5*M);
% ... non Masing Verhalten
api = para(16+5*M:16+6*M);
bp = para(17+6*M);
% ... Dehnungsgedächtnissfläche
% eta = para(18+6*M);
% cme = para(19+6*M);
% omega = para(20+6*M);
% ... proportionales Ratchetting
Qi = para(21+6*M:20+7*M);
achi = para(21+7*M);
bchi = para(22+7*M);
% ... nichtproportionales Ratchetting
cchii = para(23+7*M:22+8*M);
% ... default Korrekturterme
[aki,bk] = defaultak(M);
b = 200; % default geschwindigkeit für annäherung der radien an die Target Werte

% -------------------------------------------------------------------------
% Teilbackstress
% -------------------------------------------------------------------------
normAi = sqrt( sum( (PLINE * Ai ) .* Ai ));
normAi( normAi == 0) = 1e-40;
Li = Ai./normAi;

% -------------------------------------------------------------------------
% Berechne Exponenten chi
% -------------------------------------------------------------------------
% ... chi0
chii0 = Qi .* ( 1 + ( achi/(1+bchi*rm)^2) );
% ... chi
chii = chii0 + (chii0+0.1).*(cchii-1).*(1-abs(n' * PCHECK * Li));
chii( chii <= 0 ) = 0;

% -------------------------------------------------------------------------
% Parameterfunktion/ Target Values der Radien
% -------------------------------------------------------------------------
% ... nichtproportional
qn = qn0 * (1+an/(1+bn*rm)^2);
% ... proportional
qpi = 1+api./(1+bp*rm)^2;
% ... Endwert radius Fließfläche
rinfi_0 = rinfi(1) * ( A*qn + (1-A)*qpi(1) );
% ... Endwert radius gedächtnissfläche
rinfi_i = rinfi(2:end) .* qpi(2:end) .* ( 1 + aki./(1+bk*chii) );          % So stehts im Code
% rinfi_i = rinfi(2:end) .* qpi(2:end) .* ( 1 + aki./(1+bk*chii).^2 );          % So stehts in der Diss.
% ... Targetvalues
dummy =  1 + a1i./(1+b1*p)^2 + a2i./(1+b2*p)^2 + a3i./(1+b3*p)^2;
rT0 = rinfi_0 * dummy(1);
rTi = rinfi_i .* dummy(2:end);

% -------------------------------------------------------------------------
% Ableitungen ri i=0...M
% -------------------------------------------------------------------------
% ... Nuller abfangen
ri(ri== 0) = 1e-40;
% ... Grenzeradius FF
dr0dp = rT0 - r0; 
if abs(dr0dp) <= 1
    dr0dp = b*dr0dp;
else
    dr0dp = b*abs(dr0dp)*dr0dp;
end
% ... Grenzradien Backstresstensoren
dridp = rTi' - ri;
idx = abs(dridp) <= 1;
dridp(idx) = b*dridp(idx); 
dridp(~idx)= b*abs(dridp(~idx)).*dridp(~idx);

% -------------------------------------------------------------------------
% Ableitungen Backstresstensoren
% -------------------------------------------------------------------------
% ... Wichtungsfunktionen
Wi = normAi./abs(ri');
Wi(Wi > 1) = 1;
Wi = Wi.^chii;
% ... Schleife über Teilbackstresstensoren
dalphadp = zeros(ntens,M);
dummy = AHAT * n;
for i = 1:M 
    dalphadp(:,i) = ci(i) * ( ri(i)*dummy - Wi(i)*Ai(:,i) ) + ...
        Ai(:,i)./ri(i)*dridp(i);
end

% dummy = PHAT * n;
% dalphadp = ci.* ( ri.*dummy - Wi.*alpha) + alpha./ri.*dridp;

% ... gesamtbackstresstensor
dadp = sum(dalphadp,2);

% -------------------------------------------------------------------------
% Plastischer Tangentenmodul
% -------------------------------------------------------------------------
h = n' * PCHECK * dadp + dr0dp;
% ... begrenzungen durch E Modul
if h < -0.6 *E
    h = -0.6*E;
    dr0dp = h - nTrans * PCHECK * dadp;
end
end % Ende plastischer Tangentenmodul































% ----------------------------------------------------------------------- %
%       elastisch-plastischer tangentiale Steifigkeitstensor              %
%          (kontinuumsmechanisch Korrekt)                                 %
% ----------------------------------------------------------------------- %
function CEP = TangSteifEP(ntens,M,ink_flag,X,para,...
                           C,P,PHAT,PLINE,PCHECK,AHAT)
% Funktion berechnet den elastisch-plastischen Tangentialen
% Steifigkeitstensor 
% 
% INPUT:
%  ntens       - Spannungszustand
%  M           - Anzahl Backstress
% ink_flag     - Laststeuerung
%  X           - Zustandsvariablen
%  para        - Materialparameter
%
%  C           - elastischer Steifigkeitstensor
%  P           - Abbildung Spannung auf Spannungsdeviator
%  PHAT        - Abbildung Spannung auf Dehnung ( verdoppelt ND
%                Komponenten = PLINE für 3D)
%  PLINE       - Skalarprodukt Spannungsdeviatoren
%  PCHECK      - Skalarprodukt dev Spannung, dev Dehnung (=EINS für 3D)
%  AHAT        - Dehnunen auf Spannungen (elemenieren doppelter ND
%                Komponenten)
%
% OUTPUT:
%  CEP         - Tangentialer Steifigkeitstensor (Voigt Notation)
%__________________________________________________________________________
% -------------------------------------------------------------------------
% Identifikation der Zustandsvariablen
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteurung
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps-epsp);
else             % Dehnungssteurung
    sig = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
end
% ... Backstress
alpha = reshape(X(2*ntens+1:(2+M)*ntens),ntens,M);
a = sum(alpha,2);
% ... Backstrain
beta = X((2+M)*ntens+1:(3+M)*ntens);
% ... akk. plast. Dehn
p = X((3+M)*ntens+1);
% ... Gedäch.Fläche Dehn.raum
rm = X((3+M)*ntens+2);
% ... NP Para
A = X((3+M)*ntens+3);
% ... Radius FF in dev. Spannungsraum
r0 = X((3+M)*ntens+4);
% ... Grenzradien Backstresstensoren in dev. Spanungsraum
ri = X((3+M)*ntens+5 : (3+M)*(ntens+1) + 1);
% ... Tanaka Tensors ( als Vektor gespeichert )
% CT = X((3+M)*(ntens+1) + 2 : (3+M+ntens/2)*(ntens+1) + 1);

% -------------------------------------------------------------------------
% Fließflächennormale
% -------------------------------------------------------------------------
% ... Spannungsdeviator
s = P*sig;
% ... effektivspannung
eff = s - a;
normeff = sqrt( sum( (PLINE*eff).*eff ));
% ... Fließflächennormale (als Dehnung mit dopplten Nebendiagonalelementen)
n = PHAT* eff./normeff;
nTrans = n'; % Transposition weil die Später dauernt gebraucht wird

% -------------------------------------------------------------------------
% plastischer Tangentenmodul
% -------------------------------------------------------------------------
h = tangentenmodul(ntens,M,para,alpha,p,r0,ri,rm,A,n,...
                            PLINE,PCHECK,AHAT);
                        
% -------------------------------------------------------------------------
% Zwischenprodukte
% -------------------------------------------------------------------------
% ... C_ijkl : n_kl
Cn = C * n;
% ... n_ij : C_ijkl : n_kl
nCn = nTrans * Cn;
% ... (C_ijkl : n_kl) x  (C_mnop : n_op)
CnCn = Cn * Cn';

% -------------------------------------------------------------------------
% Tangentialer Steifigkeitstensor
% -------------------------------------------------------------------------
CEP = C - CnCn./(h + nCn);
end % Ende Tangentialer Steifigkeitsmatrix































% ----------------------------------------------------------------------- %
%       elastisch-plastischer tangentialer Nachgebigkeitstensor           %
%          (kontinuumsmechanisch Korrekt)                                 %
% ----------------------------------------------------------------------- %
function DEP = TangNachEP(ntens,M,ink_flag,X,para,...
                           D,P,PHAT,PLINE,PCHECK,AHAT)
% Funktion berechnet den elastisch-plastischen Tangentialen
% Nachgebigkeitstensor 
% 
% INPUT:
%  ntens       - Spannungszustand
%  M           - Anzahl Backstress
% ink_flag     - Laststeuerung
%  X           - Zustandsvariablen
%  para        - Materialparameter
%
%  D           - elastischer Nachgebigkeitstensor
%  P           - Abbildung Spannung auf Spannungsdeviator
%  PHAT        - Abbildung Spannung auf Dehnung ( verdoppelt ND
%                Komponenten = PLINE für 3D)
%  PLINE       - Skalarprodukt Spannungsdeviatoren
%  PCHECK      - Skalarprodukt dev Spannung, dev Dehnung (=EINS für 3D)
%  AHAT        - Dehnunen auf Spannungen (elemenieren doppelter ND
%                Komponenten)
%
% OUTPUT:
%  DEP         - Tangentialer Nachgebigkeitstensor  (Voigt Notation)
%__________________________________________________________________________
% -------------------------------------------------------------------------
% Identifikation der Zustandsvariablen
% -------------------------------------------------------------------------
if ink_flag == 0 % Spannungssteurung
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps-epsp);
else             % Dehnungssteurung
    sig = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
end
% ... Backstress
alpha = reshape(X(2*ntens+1:(2+M)*ntens),ntens,M);
a = sum(alpha,2);
% ... Backstrain
beta = X((2+M)*ntens+1:(3+M)*ntens);
% ... akk. plast. Dehn
p = X((3+M)*ntens+1);
% ... Gedäch.Fläche Dehn.raum
rm = X((3+M)*ntens+2);
% ... NP Para
A = X((3+M)*ntens+3);
% ... Radius FF in dev. Spannungsraum
r0 = X((3+M)*ntens+4);
% ... Grenzradien Backstresstensoren in dev. Spanungsraum
ri = X((3+M)*ntens+5 : (3+M)*(ntens+1) + 1);
% ... Tanaka Tensors ( als Vektor gespeichert )
% CT = X((3+M)*(ntens+1) + 2 : (3+M+ntens/2)*(ntens+1) + 1);

% -------------------------------------------------------------------------
% Fließflächennormale
% -------------------------------------------------------------------------
% ... Spannungsdeviator
s = P*sig;
% ... effektivspannung
eff = s - a;
normeff = sqrt( sum( (PLINE*eff).*eff ));
% ... Fließflächennormale (als Dehnung mit dopplten Nebendiagonalelementen)
n = PHAT* eff./normeff;
nTrans = n'; % Transposition weil die Später dauernt gebraucht wird

% -------------------------------------------------------------------------
% plastischer Tangentenmodul
% -------------------------------------------------------------------------
h = tangentenmodul(ntens,M,para,alpha,p,r0,ri,rm,A,n,...
                            PLINE,PCHECK,AHAT);
                        
% -------------------------------------------------------------------------
% tangentiale Nachgebigkeit
% -------------------------------------------------------------------------
DEP = D + (n * nTrans)/h; 
end % Ende Nachgebigkeitstensor 