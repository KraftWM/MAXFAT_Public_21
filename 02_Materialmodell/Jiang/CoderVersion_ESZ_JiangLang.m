%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung ESZ Jiang Modell für Integration nach Lang           %
%                                                                         %
%    Aufgerufen in:                                                       %
%    jianglang.m                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = CoderVersion_ESZ_JiangLang(~, X,M,eM, desig, para, epara, CEL, DEL, ...
                            P, P_line, P_hat, A, P_check)
% Inkrements des System von DGL fürs Jiang Modell mit den Ansatz von Lang
%
% INPUT:
%      ~       -> Zeit wird hier nicht gebraucht, aber wird von ode solver
%                 benötigt
%     X        -> aktueller Zustand
%    M         -> Anzahl backstresstensoren
%  para        -> Material- und Struckturparameter
%  desig       -> inkrement des pseudospannungstensors
% CEL,DEL      -> Elastische Steifigkeit und Nachgiebigkeit   
% P, P_line... -> Diverse Abbildungen
%
% OUTPUT:
%   dX   -> Inkrement der Zustandsgrößen
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%    Parameter
% -------------------------------------------------------------------------

% Spannungszustand
ntens = 3;
% Anzahl Backstresstensoren
% M = (length(para)-14)/14;                                                   
% Material
a_chi = para(3);
b_chi = para(4);
cm = para(8);
c_i0 = para(9:9+M-1);
a_i1 = para(9+M:9+2*M-1);
b_i1 = para(9+2*M:9+3*M-1);
a_i2 = para(9+3*M:9+4*M-1);
b_i2 = para(9+4*M:9+5*M-1);
r_i = para(9+5*M:9+6*M-1);
Q_i = para(9+6*M:9+7*M-1);
% Strucktur
ea_chi = epara(3);
eb_chi = epara(4);
eak = epara(5);
eck = epara(6);
ek1 = epara(7);
ecm = epara(8);
ec_i0 = epara(9:9+eM-1);
ea_i1 = epara(9+eM:9+2*eM-1);
eb_i1 = epara(9+2*eM:9+3*eM-1);
ea_i2 = epara(9+3*eM:9+4*eM-1);
eb_i2 = epara(9+4*eM:9+5*eM-1);
er_i = epara(9+5*eM:9+6*eM-1);
eQ_i = epara(9+6*eM:9+7*eM-1);

%--------------------------------------------------------------------------
%            Zustandsvariablen
%--------------------------------------------------------------------------

% Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
eeps = X(1:ntens);
% "reale" Plastische Dehnungen 
epsp = X(ntens+1:2*ntens);
% pseudo spannungen 
esig = CEL * (eeps - epsp);
% pseudo Backstress
ealphai = reshape(X(2*ntens+1 : (eM+2)*ntens),3,eM);
% "reale" plastische Bogenlänge
p = X((eM+2)*ntens+1);
% Pseudo memory surface
erm = X((eM+2)*ntens+2);
% "reale" Backstresstensoren
alphai = reshape(X((eM+2)*ntens+3 : (eM+M+2)*ntens+2),3,M);
% "reale" memory surface
rm = X((eM+M+2)*ntens+3);
% numerisch null
delta = 1e-40;
% Abfangen von Nullen 
if rm == 0
    rm = delta;
end
if erm == 0
    erm = delta;
end
% radius Pseudo Fließfläche
ek = ek1 .* ( 1 + eak * exp( eck * erm ) );

%--------------------------------------------------------------------------
%               Normale an die Strucktur Fließfläche
%--------------------------------------------------------------------------

% pseudobackstress
ea = sum(ealphai,2);
normea = sqrt( sum( (P_line * ea) .* ea ) );
if normea == 0
    normea = delta;
end

% backstress
a = sum(alphai,2);
norma = sqrt( sum( (P_line * a) .* a ) );
if norma == 0
    norma = delta;
end
% aktueller pseudo Spannungsdeviator
es = P * esig;
% pseudo effektive spannung
ebeta = es  - ea;

% normale
en = P_hat * (ebeta./(sqrt(2)*ek));
transen = en';

%--------------------------------------------------------------------------
%               Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% normen der teilbackstresstensoren
norm_ai = sqrt( sum( ( P_line * alphai ) .* alphai ) );
norm_ai(norm_ai == 0) = delta;

% normen der pseudo teilbackstresstensoren
norm_eai = sqrt( sum( ( P_line * ealphai ) .* ealphai ) );
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

% Materialparameter 
c_i = c_i0 .* (ones(1,M) + a_i1.*exp(-b_i1.*p) + a_i2.*exp(-b_i2.*p));
ec_i = ec_i0 .* (ones(1,eM) + ea_i1.*exp(-eb_i1.*p) + ea_i2.*exp(-eb_i2.*p));

chi_i = Q_i .* (2 - transen * P_check * Li) .* (1 + a_chi * exp(b_chi * rm));
echi_i = eQ_i .* (2 - transen * P_check * eLi) .* (1 + ea_chi * exp(eb_chi * erm));

% weitere hilfsgrößen
var1 = c_i .* r_i;
evar1 = ec_i .* er_i;
var2 = A * en;
var3(var3>1) = 1;
evar3(evar3>1) = 1;

% init ableitungen
dealphai_dp = zeros(ntens,eM);
dalphai_dp = zeros(ntens,M);

% Schleife über alle Backstresstensoren
% !!!! BERECHNUNGEN WÜRDE AUCH OHNE SCHLEIFE GEHEN ABER OHNO WANG HAT
% GEZEIGT; DASS ES MIT SCHLEIFE BISSI SCHNELLER IS !!!!
for ii = 1 : M
    
    % backstress
    dalphai_dp(:,ii) = var1(ii) * (var2 - (var3(ii)).^(chi_i(ii)+1).*Li(:,ii));
    
end

for ii = 1 : eM
       
    % pseudo backstress
    dealphai_dp(:,ii) = evar1(ii) * ( var2 - (evar3(ii)).^(echi_i(ii)+1).*eLi(:,ii) );
    
end

da_dp = sum(dalphai_dp,2);
dea_dp = sum(dealphai_dp,2);

%--------------------------------------------------------------------------
%                  Ableitungen der Gedächtniss- und Fließflächen                     
%--------------------------------------------------------------------------

% Hilfsgrößen
L = a./norma;
eL = ea./normea;

dummy = sign(norma - rm);
hg = 0.5 * (dummy + abs(dummy));

dummy = sign(normea - erm);
ehg = 0.5 * (dummy + abs(dummy));

dummy = sum( (P_line * da_dp) .* L);
mla = 0.5 * ( dummy + abs(dummy) );

dummy = sum( (P_line * dea_dp) .* eL);
emla = 0.5 * ( dummy + abs(dummy) );

dummy = 1 - norma/rm;
dummy = 0.5 * ( dummy + abs(dummy) );

edummy = 1 - normea/erm;
edummy = 0.5 * ( edummy + abs(edummy) );

% Ableitung Gedächtnissflächen
drm_dp = hg * mla - cm * dummy;
derm_dp = ehg * emla - ecm * edummy;

% Ableitungen des radius
dek_dp = ek1 * eak * eck * exp(eck*erm) * derm_dp;
% Abfangen von 0*Inf = NaN
dek_dp(isnan(dek_dp)) = 0; 

%--------------------------------------------------------------------------
%                  plastischer tangentenmodul                   
%--------------------------------------------------------------------------

eh = transen * P_check * dea_dp + sqrt(2) * dek_dp;

%--------------------------------------------------------------------------
%                  inkrement plastische Bogenlänge                   
%--------------------------------------------------------------------------

dp = (transen * desig)/eh;

%--------------------------------------------------------------------------
%                  inkremente der zustandsvariablen                  
%--------------------------------------------------------------------------
depsp = en .* dp;
deeps = depsp + DEL * desig;
dalphai = dalphai_dp .* dp;
drm = drm_dp .* dp;
dealphai = dealphai_dp .* dp;
derm = derm_dp .* dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente                        
%-------------------------------------------------------------------------- 
  
dX = [deeps; depsp; reshape(dealphai,3*eM,1); dp; derm; reshape(dalphai,3*M,1); drm];

end

