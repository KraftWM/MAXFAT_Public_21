function dZVAR = ...
          CoderVersion_OWT_ESZ_Lang(~, ZVAR,M, eM, desig, para, epara, CEL, DEL, ...
                            P, P_line, P_hat, A, P_check)
% Konkretes Materialmodell für ESZ und Ansatz nach Lang
%
%   INPUT:
%      ~       -> Zeit wird hier nicht gebraucht, aber wird von ode solver
%                 benötigt
%     X        -> aktueller Zustand
%  para        -> Material- und Struckturparameter
%  desig       -> inkrement des pseudospannungstensors
% CEL,DEL      -> Elastische Steifigkeit und Nachgiebigkeit   
% P, P_line... -> Diverse Abbildungen
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
% M = (length(para)-20)/6; % TODO

% Materialmodell
% kinematische Verfestigung
c_i = para(3:2+M);
r_i = para(3+M:2+2*M);
chi_i = para(3+2*M:2+3*M);
% NP Verfestigung
gamma = para(3+3*M);
gamma_np = para(4+3*M);
gamma_a = para(5+3*M);
gamma_c = para(6+3*M);
Qnpmax = para(7+3*M);
eta = para(8+3*M);
omega = para(9+3*M);
cg = para(10+3*M);

% Strukturmodell
% kinematische Verfestigung
ec_i = epara(3:2+eM);
er_i = epara(3+eM:2+2*eM);
echi_i = epara(3+2*eM:2+3*eM);
% NP Verfestigung
egamma = epara(3+3*eM);
egamma_np = epara(4+3*eM);
egamma_a = epara(5+3*eM);
egamma_c = epara(6+3*eM);
eQnpmax = epara(7+3*eM);
eeta = epara(8+3*eM);
eomega = epara(9+3*eM);
ecg = epara(10+3*eM);
% Radius FF
er0 = epara(11+3*eM);
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
delta = 1e-40;                                                             % numerisch 0

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
eeps = ZVAR(1:ntens);
% "reale" Plastische Dehnungen 
epsp = ZVAR(ntens+1:2*ntens);
% pseudo spannungen 
esig = CEL * (eeps - epsp);
% pseudo Backstress
ealphai = reshape(ZVAR(2*ntens+1 : (eM+2)*ntens),ntens,eM);
% "reale" plastische Bogenlänge
p = ZVAR((eM+2)*ntens + 1);
% Zusätzlicher Radius FF für NPV für Strukturmodell
eQ = ZVAR((eM+2)*ntens + 2);
% (Backstraintensor) Mittelpunkt Gedächtnisfläche im Dehnungsraum für Strukturmodell
ebeta = ZVAR((2+eM)*ntens+3:(3+eM)*ntens+2);
% Radius Dehnungsgedächtnisfläche für Strukturmodell
eq = ZVAR((3+eM)*ntens+3);
% Nichtproportionaalitätskennwert für Strukturmodell
eANP = ZVAR((3+eM)*ntens+4);
% Nichtproportionalitätstensor nach Tanaka für Strukturmodell
eCT = ZVAR((3+eM)*ntens+5 : (3+eM+(ntens+1)/2)*ntens+4);
% "reale" Backstresstensoren
alphai = reshape( ZVAR( (3+eM+(ntens+1)/2)*ntens+5 : (3+eM+M+(ntens+1)/2)*ntens+4 ),ntens,M);
% Zusätzlicher Radius FF für NPV 
Q = ZVAR((3+eM+M+(ntens+1)/2)*ntens+5);
% (Backstraintensor) Mittelpunkt Gedächtnisfläche im Dehnungsraum 
beta = ZVAR((3+eM+M+(ntens+1)/2)*ntens+6:(4+eM+M+(ntens+1)/2)*ntens+5);
% Radius Dehnungsgedächtnisfläche
q = ZVAR((4+eM+M+(ntens+1)/2)*ntens+6);
% Nichtproportionaalitätskennwert
ANP = ZVAR((4+eM+M+(ntens+1)/2)*ntens+7);
% Nichtproportionalitätstensor nach Tanaka
CT = ZVAR((4+eM+M+(ntens+1)/2)*ntens+8 : (4+eM+M+(ntens+1))*ntens+7);

% -------------------------------------------------------------------------
%                        Radius der Fließfläche
% -------------------------------------------------------------------------
er = er0 + eQ;


%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
es = P * esig;
% Backstress
ea = sum(ealphai,2);
% normale an Fließfläche
en = (w3d2/er).* P_hat * (es-ea);
transen = en.';


%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Normen der Teilbackstresstensoren
norm_ai = sqrt( sum( (P_line*alphai) .* alphai ) );
norm_ai(norm_ai == 0) = delta;                                             % Abfangen von 0

% normen der pseudo teilbackstresstensoren
norm_eai = sqrt( sum( (P_line * ealphai) .* ealphai ) );
norm_eai(norm_eai == 0) = delta;

% hilfsgrößen
Li = zeros(3,M);
eLi = zeros(3,eM);
var3 = zeros(1,M);
evar3 = zeros(1,eM);
for kk = 1 : M
    Li(:,kk) = alphai(:,kk)/norm_ai(kk);
    var3(kk) = norm_ai(kk)/r_i(kk);
end
for kk = 1 : eM
    eLi(:,kk) = ealphai(:,kk)/norm_eai(kk);
    evar3(kk) = norm_eai(kk)/er_i(kk);
end

% weitere hilfsgrößen ( !!!! IMPLEMENTIERUNG DER SCHLEIFE SCHEINT SCHNELLER
% ZU SEIN !!!! )
var1 = c_i .* r_i;
evar1 = ec_i .* er_i;
var2 = A * en;
var3(var3>1) = 1;
evar3(evar3>1) = 1;
var4 = transen * P_check * Li;
var4 = 0.5 * (var4 + abs(var4));
evar4 = transen * P_check * eLi;
evar4 = 0.5 * (evar4 + abs(evar4));

% Schleife über alle Backstresstensoren
% init ableitungen
dealphai_dp = zeros(ntens,eM);
dalphai_dp = zeros(ntens,M);
for ii = 1 : M
    
    % backstress
    dalphai_dp(:,ii) = var1(ii) * (var2 - (var3(ii)).^(chi_i(ii)+1).*var4(ii).*Li(:,ii));
    
end
for ii = 1 : eM
    
    % pseudo backstress
    dealphai_dp(:,ii) = evar1(ii) * ( var2 - (evar3(ii)).^(echi_i(ii)+1).*evar4(ii).*eLi(:,ii) );
    
end


% Ableitung gesamter pseudo backstress
dea_dp = sum(dealphai_dp,2);

% -------------------------------------------------------------------------
%                   Gedächtnissfläche
% -------------------------------------------------------------------------
% Effektive plastische Dehnung
effs = epsp-beta;
norm2es = effs'*(A*P_check)*effs;

% Effektive plastische Dehnung Strukturmodell
eeffs = epsp-ebeta;
norm2ees = eeffs'*(A*P_check)*eeffs;

% Gedächtnissfläche
g = norm2es - q^2;

% Gedächtnissfläche Strukturmodell
eg = norm2ees - eq^2;

% Hilfsvariable ( H(g) )
Hg = 0.5 * (sign(g) + 1);
Heg = 0.5 * (sign(eg) + 1);

% Normale an die gedächtnissfläche
if norm2es == 0
    norm2es = delta;
end

if Hg == 0 
    % ... ||n*|| < 1
%     ng = effstrain./sqrt(norm2es);
    ng = effs./q;%effstrain./sqrt(norm2es);
    dqdp = - (cg * q)^omega; 
else
    % ... ||n*|| = 1
    if norm2es > delta
        ng = effs./sqrt(norm2es);
    else
        ng = en;
    end
    dqdp = eta;
end

% Normale an die gedächtnissfläche Strukturmodell
if norm2ees == 0
    norm2ees = delta;
end

if Heg == 0 
    % ... ||n*|| < 1
%     ng = effstrain./sqrt(norm2es);
    neg = eeffs./eq;
    deqdp = - (ecg * eq)^eomega; 
else
    % ... ||n*|| = 1
    if norm2ees > delta
        neg = eeffs./sqrt(norm2ees);
    else
        neg = en;
    end
    deqdp = eeta;
end

% Hilfsvariabel <ne:n>
nn = transen * (A*P_check) * ng;
nn = 0.5 * (nn + abs(nn));
nen = transen * (A*P_check) * neg;
nen = 0.5 * (nen + abs(nen));

% Isotrope Verfestigung
dqdp = dqdp*nn;
deqdp = deqdp*nen;

% Kinematische Verfestigung
dbetadp = (1-eta) * Hg * nn * ng;
debetadp = (1-eeta) * Heg * nen * neg;

% -------------------------------------------------------------------------
%                   Nichtproportionale Verfestigung
% -------------------------------------------------------------------------
Qnpinf = Qnpmax * (1 - exp(- gamma_np * q));
Qnp = ANP*Qnpinf;
dQdp = gamma * ( Qnp - Q );

eQnpinf = eQnpmax * (1 - exp(- egamma_np * eq));
eQnp = eANP*eQnpinf;
deQdp = egamma * ( eQnp - eQ );



%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

eh = transen * P_check * dea_dp + sqrt(2/3) * deQdp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------

dp = (1/eh) * (transen * desig); 

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

% Materialmodell
depsp= dp .* en;
dalphai = dalphai_dp .* dp;
dQ = dQdp * dp;
dq = dqdp * dp;
dbeta = dbetadp * dp;
[dANP,dCT] = inkTanaka(ntens,ANP,CT,en,gamma_a,gamma_c,dp);

% Strukturmodell
deeps = depsp + DEL * desig;
dealphai = dealphai_dp .* dp;
deQ = deQdp * dp;
deq = deqdp * dp;
debeta = debetadp * dp;
[deANP,deCT] = inkTanaka(ntens,eANP,eCT,en,egamma_a,egamma_c,dp);

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------


dZVAR = [deeps; depsp; reshape(dealphai,ntens*eM,1);dp;deQ;debeta;deq;deANP;deCT;...
                       reshape(dalphai,ntens*M,1);    dQ; dbeta; dq; dANP; dCT];   

end % Ende Modell ESZ
