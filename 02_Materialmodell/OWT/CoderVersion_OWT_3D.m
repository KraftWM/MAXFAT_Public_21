function dZVAR = ...
          CoderVersion_OWT_3D(~, ZVAR,M, ink, ink_flag, parameter , C, D, ...
                              P, P_hat, P_tilde)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 3D Spannungszustand                 %
%                  OWT Modell umgeschrieben für MEX Funktione
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% M = (length(parameter)-11)/3;   % TODO                                            % Anzahl TBST

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
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
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

Li = zeros(ntens,M);
var3 = zeros(1,M);
for kk = 1 : M
    Li(:,kk) = alpha(:,kk)/norm_ai(kk);
    var3(kk) = norm_ai(kk)/r_i(kk);
end
var1 = c_i .* r_i;
var2 = P_tilde * n;
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
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
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
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
end
dZVAR = [dout; depsp ; reshape(dalpha,ntens*M,1); dp; dQ; dbeta; dq; dA; dCT];   

end % Ende Modell 3D