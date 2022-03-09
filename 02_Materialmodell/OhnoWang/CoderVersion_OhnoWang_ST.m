function dX = ...
          CoderVersion_OhnoWang_ST(~, X,M, ink, ink_flag, parameter, C, D, P, ...
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
% M = (length(parameter)-3)/3;                                               % anzahl Backstresstensoren

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
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
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
% norm_ai = sqrt( sum( (P_line.*alpha) .* alpha ) );
% norm_ai(norm_ai == 0) = delta;                                             % Abfangen von 0
% % Hilfsgrößen
Li = zeros(ntens,M);
var3 = zeros(1,M);
var4 = zeros(1,M);
for kk = 1 : M
    norm_ai = sqrt( sum((P_line.* alpha(:,kk)) .* alpha(:,kk)));
    if norm_ai == 0
        norm_ai = delta;
    end
    Li(:,kk) = alpha(:,kk)/norm_ai;
    var3(kk) = norm_ai/r_i(kk);
    var4(kk) = transn * (P_check .* Li(:,kk));
end
var1 = c_i .* r_i;
var2 = A .* n;
var3(var3>1) = 1;
% var4 = transn * (P_check .* Li);
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
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
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
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
end
% dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dp]; 
% dX = [dout; depsp ; dalpha([1:M*ntens]).'; dp]; 
dX = [dout; depsp ; dalpha(:); dp];

end % Ende Modell Sigma-tau Spannungszustand