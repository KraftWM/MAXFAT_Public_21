function [tf,gf,bg,cg] = GleitungsWL(sf,ef,b,c,hyp,varargin)
% Schätzt zyklische Kurve für Gleitungen aus zyklischer Kurve für Dehnungen
% nach verschiedenen Hypothesen.
% sf,ef,b,c  - Parameter DWL
% hyp        - Vergleichshypothese
%              1 - Mises
%              2 - Tresca
%              3 - Max. Hauptdehnung (varargin = nu)
% -------------------------------------------------------------------------

% Steigeungen
bg = b;
cg = c;

% Koef
if hyp == 1
    tf = sf/sqrt(3);
    gf = sqrt(3) * ef;
elseif hyp == 2
    tf = sf/2;
    gf = 3/2*ef;
elseif hyp == 3
    nu = varargin{1};
    tf = sf/(1+nu);
    gf = 2 * ef;
end

end