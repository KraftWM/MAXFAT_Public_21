function [Kg,ng] = GleitungRambergOsgood(Ks,ns,hyp)
% Schätzt zyklische Kurve für Gleitungen aus zyklischer Kurve fuer Dehnungen
% nach verschiedenen Hypothesen.
% Ks,ns      - Parameter Ramberg Osgood für Dehnungen
% hyp        - Vergleichshypothese
%            1 - Mises
%            2 - Tresca
%            3 - Max. Hauptdehnung (varargin = nu)
% -------------------------------------------------------------------------

% Exponent 
ng = ns;

% Steifigkeit
if hyp == 1
    Kg = Ks * (1/sqrt(3))^(ns+1);
elseif hyp == 2
    Kg = 0.5 * Ks * (2/3)^ns;
elseif hyp == 3
    msg = 'net implementiert';
    error(msg);
end

end