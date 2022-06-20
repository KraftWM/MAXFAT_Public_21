%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung ESZ Ohno Wang Modell für Integration mit             %
%    Energieinkrement (nach mir)                                          %
%                                                                         %
%    Aufgerufen in:                                                       %
%    ohnowangenergie.m                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = CoderVersion_ESZ_OhnoWangEnergie(~, ZVAR, M, domega, para, ...
                       GINV, PHAT,...
                       M_, MLINE, MHAT, A, MCHECK)
% Inkrement des OhnoWang Modells bei energiegesteuerter Integration                   
%
% INPUT:
%      ~       -> Zeit wird hier nicht gebraucht, aber wird von ode solver
%                 benötigt
%  ZVAR        -> aktueller Zustand
%   M          -> Anzahl der Backstresstensoren (auf 20 beschränkt)
%  para        -> Material- und Struckturparameter
%  domega      -> inkrement der energien
%   GINV       -> inverse Ableitung Energie nach Spannungen
%  PHAT        -> ABleitung Energie nach plastischen Dehnungen (!als vektor)
% M_, MLINE... -> Diverse Abbildungen
%
% OUTPUT:
%   dX   -> Inkrement der Zustandsgrößen                   
%
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
% M = (length(para)-3)/3;                                               % anzahl Backstresstensoren

r0 = para(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = para(3:2+M);
r_i = para(3+M:2+2*M);
chi_i = para(3+2*M:2+3*M);

% oft verwendete Konstanden
w3d2 = sqrt(3/2);

delta = 1e-40;                                                             % numerisch 0

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Spannungen
sig = ZVAR(1:ntens);
% backstresstensoren
alpha = reshape(ZVAR(2*ntens+1:(M+2)*ntens),ntens,M);

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = M_ * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r0).* MHAT * (s-a);
transn = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Normen der Teilbackstresstensoren
norm_ai = sqrt( sum( (MLINE*alpha) .* alpha ) );
norm_ai(norm_ai == 0) = delta;                                             % Abfangen von 0

% Hilfsvariablen 
Li = zeros(3,M);
var3 = zeros(1,M);
for kk = 1 : M
    Li(:,kk) = alpha(:,kk)/norm_ai(kk);
    var3(kk) = norm_ai(kk)/r_i(kk);
end
var1 = c_i .* r_i;
var2 = A * n;
var3(var3>1) = 1;
var4 = transn * MCHECK * Li;
var4 = 0.5 * (var4 + abs(var4));
% Schleife über alle Backstresstensoren (aus irgend einem Grund ist
% Schleife schneller als vektorisiert)
% init ableitungen
dalpha_dp = zeros(ntens,M);
for ii = 1 : M
    
    % backstress
    dalpha_dp(:,ii) = var1(ii) * (var2 - (var3(ii)).^(chi_i(ii)+1).*var4(ii).*Li(:,ii));
    
end

% Ableitung gesamt backstress
da_dp = sum(dalpha_dp,2);

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = transn * MCHECK * da_dp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------
dummy1 = GINV *domega;
if size(PHAT,2) == 1
    dummy2 = GINV * ( PHAT .* n );
else
    dummy2 = GINV * ( PHAT * n );
end
dp = (transn * dummy1)/(h + transn * dummy2 );

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp=dp.*n;
dalpha=dalpha_dp.*dp;
dsig = dummy1 - dummy2.*dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

% dX = [dsig; depsp ; reshape(dalpha,ntens*M,1); dp]; 
dX = [dsig; depsp ; dalpha(:); dp]; 

end % Ende Modellgleichung ESZ