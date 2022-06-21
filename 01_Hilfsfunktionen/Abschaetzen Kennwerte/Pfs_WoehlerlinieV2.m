function [Pfs_WS_stuetz,Pfs_WSD_stuetz,d1,d2,kfs,f0025] = Pfs_WoehlerlinieV2(Rm,wsgruppe)
% Bestimmt die Stuetzstelle fuer Pfs, Indem in die Formel für Pfs die
% abschätzung der DehnungsWL nach der FKM Methode (Wächter) eingesetzt
% wird
%
% Pfs_WS_stuetz bei N = 100
% wsgruppe - Wrkstoffgrupe 1 = Stahl
%                          2 = Stahlguss
%                          3 = Aluknet
%                          4 = Aluguss
%                          5 = GJS
%                          6 = GJM
%                          7 = GJL
%                          8 = Hoechstfester Stahl
% Rm       - Zugefestigkeit
% -------------------------------------------------------------------------

% ... statische Werte abschätzen
[E,nu] = StatischFKM(wsgruppe);

% ... PWL
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    % ... Anteil aus SpannungsWL
    as = 3.1148;
    ak = 0.4625;
    bk = -0.897;
    bs = 0.897;
    b = -0.097;
    a1_1 = (1+nu)/E * ((2*1)^b * as + (2*1)^(2*b) * 0.5 *as^2 *ak);
    a1_100 = (1+nu)/E * ((2*100)^b * as + (2*100)^(2*b) * 0.5 *as^2 *ak);
    a1_1E5 = (1+nu)/E * ((2*1e5)^b * as + (2*1e5)^(2*b) * 0.5 *as^2 *ak);
    a1_D = (1+nu)/E * ((2*1e6)^b * as + (2*1e6)^(2*b) * 0.5 *as^2 *ak);
    % ... Anteil aus plast. DehnWL
    ef = 1033 * Rm^(-1.235);
    c = -0.52;
    if ef > 0.338
        ae = 0.338;
        be = 0;
    else
        ae = 1033;
        be = -1.235;
    end
    a2_1 = (2*1)^c * 1.5 * ae+ (2*1)^(b+c) * 0.75*ak*as*ae;
    a2_100 = (2*100)^c * 1.5 * ae+ (2*100)^(b+c) * 0.75*ak*as*ae;
    a2_1E5 = (2*1e5)^c * 1.5 * ae+ (2*1e5)^(b+c) * 0.75*ak*as*ae;
    a2_D = (2*1e6)^c * 1.5 * ae+ (2*1e6)^(b+c) * 0.75*ak*as*ae;

    f0025 = 0;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    % ... Anteil aus SpannungsWL
    as = 1.732;
    ak = 0.9006;
    bk = -0.982;
    bs = 0.982;
    b = -0.102;
    a1_1 = (1+nu)/E * ((2*1)^b * as + (2*1)^(2*b) * 0.5 *as^2 *ak);
    a1_100 = (1+nu)/E * ((2*100)^b * as + (2*100)^(2*b) * 0.5 *as^2 *ak);
    a1_1E5 = (1+nu)/E * ((2*1e5)^b * as + (2*1e5)^(2*b) * 0.5 *as^2 *ak);
    a1_D = (1+nu)/E * ((2*1e6)^b * as + (2*1e6)^(2*b) * 0.5 *as^2 *ak);
    % ... Anteil aus plast. DehnWL
    c = -0.58;
    ae = 0.847;
    be = -0.181;
    a2_1 = (2*1)^c * 1.5 * ae+ (2*1)^(b+c) * 0.75*ak*as*ae;
    a2_100 = (2*100)^c * 1.5 * ae+ (2*100)^(b+c) * 0.75*ak*as*ae;
    a2_1E5 = (2*1e5)^c * 1.5 * ae+ (2*1e5)^(b+c) * 0.75*ak*as*ae;
    a2_D = (2*1e6)^c * 1.5 * ae+ (2*1e6)^(b+c) * 0.75*ak*as*ae;

    f0025 = 0;
    
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    % ... Anteil aus SpannungsWL
    as = 9.21;
    ak = 0.1823;
    bk = -0.742;
    bs = 0.742;
    b = -0.106;
    a1_1 = (1+nu)/E * ((2*1)^b * as + (2*1)^(2*b) * 0.5 *as^2 *ak);
    a1_100 = (1+nu)/E * ((2*100)^b * as + (2*100)^(2*b) * 0.5 *as^2 *ak);
    a1_1E5 = (1+nu)/E * ((2*1e5)^b * as + (2*1e5)^(2*b) * 0.5 *as^2 *ak);
    a1_D = (1+nu)/E * ((2*1e6)^b * as + (2*1e6)^(2*b) * 0.5 *as^2 *ak);
    % ... Anteil aus plast. DehnWL
    c = -0.83;
    ae = 895.9;
    be = -1.183;
    a2_1 = (2*1)^c * 1.5 * ae+ (2*1)^(b+c) * 0.75*ak*as*ae;
    a2_100 = (2*100)^c * 1.5 * ae+ (2*100)^(b+c) * 0.75*ak*as*ae;
    a2_1E5 = (2*1e5)^c * 1.5 * ae+ (2*1e5)^(b+c) * 0.75*ak*as*ae;
    a2_D = (2*1e6)^c * 1.5 * ae+ (2*1e6)^(b+c) * 0.75*ak*as*ae;

    f0025 = 0;
else
    msg = 'Falsche Werstoffgrupppe uebergeben, Werte für Stahl werden verwendet';
    warning(msg)  
    % ... Anteil aus SpannungsWL
    as = 3.1148;
    ak = 0.4625;
    bk = -0.897;
    bs = 0.897;
    b = -0.097;
    a1_1 = (1+nu)/E * ((2*1)^b * as + (2*1)^(2*b) * 0.5 *as^2 *ak);
    a1_100 = (1+nu)/E * ((2*100)^b * as + (2*100)^(2*b) * 0.5 *as^2 *ak);
    a1_1E5 = (1+nu)/E * ((2*1e5)^b * as + (2*1e5)^(2*b) * 0.5 *as^2 *ak);
    a1_D = (1+nu)/E * ((2*1e6)^b * as + (2*1e6)^(2*b) * 0.5 *as^2 *ak);
    % ... Anteil aus plast. DehnWL
    ef = 1033 * Rm^(-1.235);
    c = -0.52;
    if ef > 0.338
        ae = 0.338;
        be = 0;
    else
        ae = 1033;
        be = -1.235;
    end
    a2_1 = (2*1)^c * 1.5 * ae+ (2*1)^(b+c) * 0.75*ak*as*ae;
    a2_100 = (2*100)^c * 1.5 * ae+ (2*100)^(b+c) * 0.75*ak*as*ae;
    a2_1E5 = (2*1e5)^c * 1.5 * ae+ (2*1e5)^(b+c) * 0.75*ak*as*ae;
    a2_D = (2*1e6)^c * 1.5 * ae+ (2*1e6)^(b+c) * 0.75*ak*as*ae;

    f0025 = 0;
    
end 

% ... berechne Stuetzstellen
kfs = ak * Rm ^bk;
Pfs_WS1_stuetz = a1_1 * Rm^bs + a2_1 * Rm^be;
Pfs_WS_stuetz = a1_100 * Rm^bs + a2_100 * Rm^be;
Pfs_WS1E5_stuetz = a1_1E5 * Rm^bs + a2_1E5 * Rm^be;
Pfs_WSD_stuetz = a1_D * Rm^bs + a2_D * Rm^be;
d1 = (log10(Pfs_WS_stuetz)-log10(Pfs_WS1_stuetz))/2;
d2 = (log10(Pfs_WS1E5_stuetz)-log10(Pfs_WS_stuetz))/(5-2);

end