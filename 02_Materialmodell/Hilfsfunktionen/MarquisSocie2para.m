function [Qnpmax,gamma_np] = MarquisSocie2para(Kstrich,nstrich,alpha)
% Bestimme Materialparameter für nichtproportionale Verfestigung für OWT
% Modell
%
% INPUT:
%     Kstrich,nstrich - Parameter Ramberg Osgood
%     alpha           - Parameter Marquis Socie Ansatz
%
% OUTPUT:
% Onpmax,gamma_np     - Parameter Plastizitätsmodell
% 
% ANMERKUNGEN:
%
% NP Fließkurve nach MS Ansatz: 
%     (1) Knp = (1+alpha/sqrt(2)) * Kstrich   
%                      1/sqrt(2) =  Tanaka Parameter für 90°
%                      Phasenverschiebung & sig_a/tau_a = sqrt(3)
%
%     (2) sigmanp = f(epsp) = Knp * epsp^nstrich     -> Fließkurve NP 
%     (3) sigma   = f(epsp) = Kstrich * epsp^nstrich -> Fließkurve Prop/Uniax
%
% Zusätzliche NP Verfestigung im plastizitätsmodell
%
%     (4) sigmanp-sigma = alpha/sqrt(2) * Kstrich * epsp^nstrich
%                       = Qnpmax(1-exp(-gamma_np*q))
%
%     (5) q = sqrt(3/2) * epsp     (ungefähr)
% -------------------------------------------------------------------------


% Bereich in dem gefittet wird
x = 0.001:0.001:0.02;          % = epsp
% Zu fittende Werte
y = alpha/sqrt(2) * Kstrich * x.^nstrich;
% Fit Function
% !!!!!!! Fit ausgestellt, Paramater werden einfach geraten
Qnpmax = alpha/sqrt(2) * Kstrich * 0.03.^nstrich;
gamma_np = 500;

end