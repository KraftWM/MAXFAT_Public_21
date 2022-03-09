function ZVAR0 = initZVARDoring(ntens,para,M)
% Initialisieren aller Zustandsvariablen des Döring Modells je nach
% Spannungsszustand
% 
%
% INPUT:
% ntens         - Anzahl Tensorkomponenten 
% para          - parameter des Döring Modells
% M             - Anzahl der Backstresstensoren
% 
% OUTPUT:
% ZVAR0         - Startwerte der Zustandsvariablen (im wesentlichen nur ri
%                 Werte)
%__________________________________________________________________________

% Anzahl Zustandsvariablen
nZV = (3+M+ntens/2)*(ntens+1)+1;

% init ZV mit nullen
ZVAR0 = zeros(nZV,1);

% auslesen parameter
Qi = para(21+6*M:20+7*M); % i = 1...M
achi = para(21+7*M);
api = para(16+5*M:16+6*M); % i = 0...M
a1i = para(4+2*M:4+3*M); % i = 0...M   
a2i = para(5+3*M:5+4*M); % i = 0...M
a3i = para(6+4*M:6+5*M); % i = 0...M
b1 = para(7+5*M);
b2 = para(8+5*M);
b3 = para(9+5*M);
rinf0i = para(3+M:3+2*M); % i = 0...M
[aki,bk] = defaultak(M);  % i = 1...M

% Setzte deviatorische Fließspannung ( r(i=0) )
sf = para(10+5*M);     % Startwert r00
r0 = rinf0i(1);        % Unendlich Wert r0
if sf < 2e-10
    p = -sf;
    sf = r0 * sqrt(3/2) * ( 1 + api(1) );
    sf = sf * ( 1 + a1i(1) / ( 1 + b1*p )^2 ...
                  + a2i(1) / ( 1 + b2*p )^2 ...
                  + a3i(1) / ( 1 + b3*p )^2 );
else
    sf = sf * sqrt(3/2);
    p = 0;
end
ZVAR0((3+M)*ntens+4) = sf;
ZVAR0((3+M)*ntens+1)= p; 

% Setze begrenzungsradien der Backstresstensoren (= Targetvalue von p=0)
% r(i=1..M)
ri = zeros(M,1);
for i = 1:M
    chi0 = Qi(i) * (1+achi);
%     dummy = aki(i)/(1+bk*chi0)^2+1;             % SO stehts in der Diss
    dummy = aki(i)/(1+bk*chi0)+1;                % So stehts im Code
    rinfi = rinf0i(i+1)*dummy*(1+api(i+1));
    ri(i) = rinfi * ( 1 + a1i(i+1) / ( 1 + b1*p )^2 ...
                        + a2i(i+1) / ( 1 + b2*p )^2 ...
                        + a3i(i+1) / ( 1 + b3*p )^2 ); 
end
ZVAR0((3+M)*ntens+5 : (3+M)*(ntens+1) + 1) = ri;
