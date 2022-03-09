function dX = CoderVersion_Jiang_ST(~, X,M, ink, ink_flag, parameter, C, D, ...
                                 P, P_line, P_hat, A, P_check)
% Konkretes Materialmodell nach Jiang für sigma-tau Spannungszustand
%
% INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2021 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------

%-----------------------------------------------------------------------
%            Identifikation Materialparameter
%-----------------------------------------------------------------------
ntens = 2;                                                             % Tensorkomponenten
% M = (length(parameter)-8)/7;                                           % Anzahl Backstresstensoren

a_chi = parameter(3);
b_chi = parameter(4);
ak = parameter(5);
ck = parameter(6);
k1 = parameter(7);
cm = parameter(8);
c_i0 = parameter(9:9+M-1);
a_i1 = parameter(9+M:9+2*M-1);
b_i1 = parameter(9+2*M:9+3*M-1);
a_i2 = parameter(9+3*M:9+4*M-1);
b_i2 = parameter(9+4*M:9+5*M-1);
r_i = parameter(9+5*M:9+6*M-1);
Q_i = parameter(9+6*M:9+7*M-1);

%-----------------------------------------------------------------------
%            Identifikation der Zustandsvariablen
%-----------------------------------------------------------------------

if ink_flag == 0
    eps = X(1:ntens);                                               % Dehnung aus Zustandsvariable lesen
    epsp = X(ntens+1:2*ntens);                                      % Plastischer Anteil der Dehnung aus Zustandsvariable lesen
    sig = C * (eps - epsp);                                         % Spannung berechnen
    
elseif ink_flag == 1
    sig = X(1:ntens);                                               % Spannung aus Zustandsvariable lesen
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')    
end
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);                              % akt. Backstresstensoren ermitteln

%-----------------------------------------------------------------------
%               Radius der Fließfläche berechnen
%-----------------------------------------------------------------------
delta = 1e-40;
p = X((M+2)*ntens+1);                                                   % plastische Bogenlänge
rm = X((M+2)*ntens+2);                                                  % Radius der Gedächtnisfläche
if rm == 0
    rm = delta;                                                        % Durch 0 teilen abfangen
end

%-----------------------------------------------------------------------
%               Normale an die Fließfläche berechnen
%-----------------------------------------------------------------------
s = P .* sig;                                                            % Deviator des Spannungstensors
a = sum(alpha,2);                                                       % Backstress
k = k1 .* ( 1 + ak * exp( ck * rm ) );                                  % radius FF
beta = s - a;                                                           % effektive Spannung

n = P_hat .* (beta./(sqrt(2)*k));
nTrans = n';                                                            % Transposition

norm_a = sqrt( sum( (P_line.*a) .* a) );                                 % Norm Backstress
norm_a( norm_a == 0) = delta;

%-----------------------------------------------------------------------
%              Ableitungen der Teilbackstresstensoren
%-----------------------------------------------------------------------

c_i = c_i0 .* (ones(1,M) + a_i1.*exp(-b_i1.*p) ...                      % Variable für zyklische Ver- bzw. Entfestigung
    + a_i2.*exp(-b_i2.*p));

% Norm von alpha_i ermitteln
% norm_ai = sqrt( sum( (P_line.*alpha) .* alpha ) );
% norm_ai(norm_ai == 0) = delta;

%    % Hilfsvariablen
Li = zeros(ntens,M);
var1 = zeros(1,M);
var0 = zeros(1,M);
for kk = 1 : M
    norm_ai = sqrt( sum((P_line.* alpha(:,kk)) .* alpha(:,kk)));
    if norm_ai == 0
        norm_ai = delta;
    end
    Li(:,kk) = alpha(:,kk)/norm_ai;
    var1(kk) = norm_ai./r_i(kk);
    var0(kk) = nTrans * (P_check .* Li(:,kk)); 
end
% Hilfsvariable
chi = Q_i .* (2 - var0) .* (1 + a_chi * exp(b_chi * rm));
var0 = c_i .* r_i;
var1( var1 > 1) = 1;
var2 = A .* n;

% schleife über alle Backstresstensoren
dalpha_dp = zeros(ntens,M);
for i = 1 : M
    
    dalpha_dp(:,i) = var0(i) * ( var2 - var1(i).^(chi(i)+1).* Li(:,i) );
    
end

% Ableitung gesamtbackstress
da_dp = sum(dalpha_dp,2);                                               % Ableitung: Backstresstensor nach plast. Bogenlänge


%-----------------------------------------------------------------------
%                  Ableitungen der Gedächtnissfläche
%-----------------------------------------------------------------------

L = a./norm_a;                                                         % Normierter Backstresstensor
dummy = sign(norm_a - rm);                                             % Hilfsvariable
hg = 0.5 * ( dummy + abs(dummy));                                      % Heaviside
dummy = sum( (P_line .* da_dp) .* L);                                   % Hilfsvariable
mla = 0.5 * ( dummy + abs(dummy) );                                    % maccaulay
dummy = 1- norm_a/rm;                                                  % Hilfsvariable
dummy = 0.5 * ( dummy + abs(dummy) );                                  % Hilfsvariable

% Ableitung Mem. surf
drm_dp = hg * mla - cm * dummy;

% Ableitung FF
dk_dp = k1 * ak * ck * exp(ck*rm) * drm_dp;                            % Ableitung: Fließfläche nach plast. Bogenlänge
% Abfangen von 0*Inf = NaN
dk_dp(isnan(dk_dp)) = 0;

%-----------------------------------------------------------------------
%                       Plastischer Modul
%-----------------------------------------------------------------------

h = nTrans * (P_check .* da_dp) + sqrt(2) * dk_dp;                        % Plastischer Modul

%-----------------------------------------------------------------------
%                Inkrement der plastischen Bogenlänge
%-----------------------------------------------------------------------

if ink_flag == 0
    
    dp = (nTrans * ink)/h;                                              % Plast. Inkrement bei Spannungssteuerung
    
elseif ink_flag == 1
    
    dp = ( nTrans * (C * ink)) / (h + (nTrans * (C * n)));              % Plast. Inkrement bei Dehnungssteuerung
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')    
end

%-----------------------------------------------------------------------
%                  Inkremente der Zustandsvariablen
%-----------------------------------------------------------------------

depsp = n * dp;                                                        % Inkrement der plast. Dehnung
dalpha = dalpha_dp * dp;                                               % Inkrement der Backstresstensoren
drm = drm_dp * dp;                                                     % Inkrement Radius der Gedächtnisfläche

%-----------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%-----------------------------------------------------------------------

if ink_flag == 0
    deps = D * ink + depsp;                                              % Inkrement der Dehnung berechnen
    out = vertcat(deps, depsp);                                        % Hilfsvariable
    
elseif ink_flag == 1
    dsig = C * (ink - depsp);                                          % Inkrement der Spannung berechnen
    out = vertcat(dsig, depsp);                                        % Hilfsvariable
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')    
end
dX = vertcat(out, reshape(dalpha,ntens*M,1), dp, drm);                                       % Inkremente der Zustandsvariablen

end % Ende sigma-tau