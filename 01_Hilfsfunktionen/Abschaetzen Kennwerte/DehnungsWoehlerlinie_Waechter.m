function [sf,ef,b,c] = DehnungsWoehlerlinie_Waechter(wsgruppe,Rm)
% Bestimmt die Parameter der Dehnungswöhlerlinie nach Manson-Coffin-Morrow
% Gleichung 
%
% Abschätzung nach Diss. Wächter Tabelle 7.28
% 
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

% Besimme Werkstoffabhänige abschätzung
if any(wsgruppe == 1) || strcmp(wsgruppe,'Stahl')
    Rmmin = 121;
    Rmmax = 2296;
    fak_sf = 3.1148;
    exp_sf = 0.897;
    fak_ef = 1033;
    exp_ef = -1.235;
    b = -0.097;
    c = -0.52;
    epsgrenz = 0.338;
elseif any(wsgruppe == 2) || strcmp(wsgruppe,'Stahlguss')
    Rmmin = 496;
    Rmmax = 1144;
    fak_sf = 1.732;
    exp_sf = 0.982;
    fak_ef = 0.847;
    exp_ef = -0.181;
    b = -0.102;
    c = -0.58;
    epsgrenz = 100;
elseif any(wsgruppe == 3) || strcmp(wsgruppe,'Aluknet')
    Rmmin = 216;
    Rmmax = 649;
    fak_sf = 9.12;
    exp_sf = 0.742;
    fak_ef = 895.9;
    exp_ef = -1.183;
    b = -0.106;
    c = -0.83;
    epsgrenz = 100;
else
    msg = ['Falsche Werstoffgrupppe uebergeben, Dehnungswöhlerlinie',...
       'kann nicht abgeschätzt werden. Werte für Stahl werden verwendet.'];
    warning(msg)    
    Rmmin = 121;
    Rmmax = 2296;
    fak_sf = 3.1148;
    exp_sf = 0.897;
    fak_ef = 1033;
    exp_ef = -1.235;
    b = -0.097;
    c = -0.52;
    epsgrenz = 0.338;
end 

% Checke Grenzen
if Rm > Rmmax || Rm < Rmmin
    msg = 'Zugfestigkeit außerhalb der Grenzen der Abschätzmethode';
    warning(msg)  
end

% Abschätzen
sf = fak_sf * Rm^exp_sf;
ef = min([epsgrenz fak_ef*Rm^exp_ef]);
end