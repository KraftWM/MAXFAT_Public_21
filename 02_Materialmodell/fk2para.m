function para = fk2para(sig, epsp, E, nu, r0, chi_i, M, modell)
% Funktion zum bestimmen der Materialparameter aus Materialfließkurve
%
% Stützstellen in plastischen Dehnungen und Spannungen in dieser Funktion
% müssen schon aufbereitet sein
%
% INPUT:
% sig    -> Stützstellen Spannungen
% epsp   -> Stützstellen Dehnungen
% E,nu   -> elastizitätskonstanten
%  r0    -> Fließspannung
% chi_i  -> Rat. para
%  M     -> Anzahl Backstress
% modell -> Materialmodell
%
%
% OUTPUT:
% para   -> Materialparameter
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Tangentenmoduli
ii = 2:length(sig);
iim1 = 1:length(sig)-1;
H = zeros(1,M+1);
H(1:M) = (sig(ii)-sig(iim1))./(epsp(ii)-epsp(iim1));


% -------------------------------------------------------------------------
% c_i und r_i
ii = 1:length(sig)-1;
iip1 = 2:length(sig);
c_i = sqrt(2/3) ./ epsp(2:end);
r_i = 2/3 .* (H(ii)-H(iip1))./c_i; r_i(r_i < 0) = 1e-40;

% -------------------------------------------------------------------------
% zuweisen der Parameter je nach Modell
switch modell
    case 'Chaboche'
        para = NaN(1,2*M+5);
        c_i = c_i * sqrt(3/2);                                                   % Anpassen an meine Implementierung
        r_i = r_i * sqrt(3/2);                                             % Anpassen an meine Implementierung
%         r_i = r_i .* [1+0.3/M:0.3/M:1.3];                                  % Pragmatische Korrektur
%         r_i(end) = r_i(end) * sqrt(3/2);                                 % Pragmatische Korrektur
        para(1) = E;                                                       % E-Modul
        para(2) = nu;                                                      % Querdehnzahl
        para(3) = 0;                                                       % q für isotrope Verfestigung
        para(4) = 0;                                                       % gamma für isotrope Verfestigung
        para(5:4+M) = c_i;                                                 % Parameter c_i
        para(M+5:2*M+4) = r_i;                                             % Parameter r_i
%         para(M*2+4) = r_i(M)*3/2;
        para(2*M+5) = r0;                                                  % Startradius FF
    case 'OhnoWang'
        para = NaN(1,3*M+3);
        para(1) = E;                                                       % E-Modul
        para(2) = nu;                                                      % Querdehnzahl
        para(3:2+M) = c_i;                                                 % Parameter c_i
        para(M+3:2*M+2) = r_i;                                             % Parameter r_i
        para(M*2+3:M*3+2) = chi_i;                                         % Rattcheting
        para(3*M+3) = r0;                                                  % Startradius FF
    case 'KarimOhno'
        para = NaN(1,3*M+3);
        c_i = c_i * sqrt(3/2);                                             % Anpassen an meine Implementierung
        r_i = r_i * sqrt(3/2);                                             % Anpassen an meine Implementierung
%         r_i(end) = r_i(end) * sqrt(3/2);                                   % Pragmatische Korrektur
        para(1) = E;                                                       % E-Modul
        para(2) = nu;                                                      % Querdehnzahl
        para(3:2+M) = c_i;                                                 % Parameter c_i
        para(M+3:2*M+2) = r_i;                                             % Parameter r_i
        if chi_i(1) == 5
            para(2*M+3:3*M+2) = 0.01 .* ones(1,M);
        else
            para(2*M+3:3*M+2) = 0.005 .* ones(1,M);
        end
        para(3*M+3) = r0;    
    case 'Jiang'
        para = NaN(1,8+7*M);
        para(1) = E;                                                       % E-Modul
        para(2) = nu;                                                      % Querdehnzahl
        para(3) = 0;                                                       % a_chi Jiang Modell
        para(4) = 1000 * 0;                                                % b_chi Jiang Modell
        para(5) = 0;                                                       % a_k   Jiang Modell
        para(6) = 500 * 0;                                                 % b_k   Jiang Modell
        para(7) = sig(1)/sqrt(3);                                          % k1    Jiang Modell
        para(8) = 100;                                                     % cm    Jiang Modell
        para(9:M+8) = c_i;                                                 % c0_i  Jiang Modell
        para(M+9:2*M+8) = zeros(1,M);                                      % a1_i  Jiang Modell         % a1_i  Jiang Modell
        para(2*M+9:3*M+8) = 100 * zeros(1,M);                              % b1_i  Jiang Modell
        para(3*M+9:4*M+8) = zeros(1,M);                                    % a2_i  Jiang Modell
        para(4*M+9:5*M+8) = 1000 * zeros(1,M);                             % b2_i  Jiang Modell
        para(5*M+9:6*M+8) = r_i;                                           % r_i   Jiang Modell
        para(6*M+9:7*M+8) = chi_i;                                         % Q_i   Jiang Modell
    case 'Doring'
        para = NaN(1,22+8*M);
        % ... Elast.
        para(1) = E;                                                       % E-Modul
        para(2) = nu;                                                      % Querdehnzahl
        % ... zyklisch Stabile Kurve
        para(3:2+M) = c_i;                                                 % ci Döring Modell (i = 1...M)
        para(3+M) = r0 * sqrt(2/3);                                        % r_inf_i_0 Döring Modell (i=0) Start/Endradius der FF in pi Ebene 
        para(4+M:3+2*M) = r_i;                                             % r_inf_i_0 Döring Modell (i=1...M)
        % ... transiente zyklische Ver und Entfestigung (default Werte)      
        para(4+2*M:4+3*M) = zeros(1,M+1);                                  % a1_i Döring Modell (i=0...M)
        para(5+3*M:5+4*M) = zeros(1,M+1);                                  % a2_i Döring Modell (i=0...M)
        para(6+4*M:6+5*M) = zeros(1,M+1);                                  % a3_i Döring Modell (i=0...M)
        para(7+5*M) = 100;                                                 % b1 Döring Modell  
        para(8+5*M) = 10;                                                  % b2 Döring Modell
        para(9+5*M) = 1;                                                   % b3 Döring Modell
        para(10+5*M) = para(3+M);                                          % r00 Döring Modell Startradius der FF in pi Ebene oder -p0 für trick mit negativem Radius                                        
        % ... nichtproportionale Zusatzverfestigung
        para(11+5*M) = 1;                                                  % qn0 Döring Modell
        para(12+5*M) = 0;                                                  % an Döring Modell
        para(13+5*M) = 1000;                                               % bn Döring Modell
        para(14+5*M) = 50;                                                 % ct Döring Modell
        para(15+5*M) = 100;                                                % ca Döring Modell
        % ... Non Masing Verhalten
        para(16+5*M:16+6*M) = zeros(1,M+1);                                % api Döring Modell (i=0...M)
        para(17+6*M) = 500;                                                % bp Döring Modell
        % ... Verhalten der Dehnungsgedächtnissfläche
        para(18+6*M) = 0.1;                                                % eta Döring Modell
        para(19+6*M) = 5;                                                  % cme Döring Modell
        para(20+6*M) = 1;                                                  % omega Döring Modell
        % ... proportionales Ratchetting
        para(21+6*M:20+7*M) = chi_i;                                       % Q_i Döring Modell (i=1...M)
        para(21+7*M) = 0;                                                  % a_chi Döring Modell
        para(22+7*M) = 500;                                                % b_chi Döring Modell
        % ... nichtproportionales Ratchetting
        para(23+7*M:22+8*M) = ones(1,M);                                   % c_chi Döring Modell (i=1...M)
    case 'OWT'
        para = NaN(1,3*M+11);
        % Statisch
        para(1) = E;                                                       % E-Modul
        para(2) = nu;                                                      % Querdehnzahl
        % Backstresstensor
        para(3:2+M) = c_i;                                                 % Parameter c_i
        para(M+3:2*M+2) = r_i;                                             % Parameter r_i
        para(M*2+3:M*3+2) = chi_i;                                         % Rattcheting
        % NP Verfestigung
        para(M*3+3) = 100;                                                 % gamma    (Standartwert)
        para(M*3+4) = 0;                                                   % gamma_np
        para(M*3+5) = 80;                                                  % gamma_a  (Standartwert)
        para(M*3+6) = 50;                                                  % gamma_c  (Standartwert)
        para(M*3+7) = 0;                                                   % Qnpmax
        para(M*3+8) = 0.1;                                                 % eta      (Standartwert)
        para(M*3+9) = 1;                                                   % omega    (Standartwert)
        para(M*3+10) = 5;                                                  % cg       (Standartwert)
        % Radius FF (Startwert)
        para(3*M+11) = r0;                                                 % Startradius FF
    otherwise
        msg = 'Materialmodell nicht impelementiert';
        error(msg);
end
end