function [EXTI,EXTII,EXTIII] = ExtremaSchnittdehnungen(EPS)
% Gibts die Extrema der Schnittdehnungen für gegebene Dehnungstensor im 
% ebenen Spannungszustand zurück.
% 
% INPUT:
%  EPS   - Dehnungstensor entweder (4x1) oder (1x4) 
% 
% OUTPUT:
%  EXTI   - Extrema in Normaldehnung
%          1. Zeile phi
%          2. Zeile psi
%          3. Zeile Min(-1) Max(1) Sattelpunkt(0)
%  EXTII  - Extrema in gamma_xz
%  EXTIII - Extrema in gamma_xy
% -------------------------------------------------------------------------

%% Input verwalten
nrows = size(EPS,1);
ncols = size(EPS,2);
if nrows + ncols ~= 5
    warning('Falscher Input in ExtremaSchnittdehnungen')
    EXTI = NaN;
    EXTII = NaN; 
    EXTIII = NaN;
    return;
end
exx = EPS(1);
eyy = EPS(2);
ezz = EPS(3);
gxy = EPS(4);
%% Normaldehnungen (ModeI)
phi = 0.5*atan(gxy/(exx-eyy));
phiI = [];
% finde alle Phi in [0 2pi]
while phi <= 2 * pi
    phiI = [phiI,phi];
    phi = phi + pi/2;
end
phiI = phiI(phiI > 0);
phiI = phiI(phiI <= 2* pi);
% Unterscheide Maxima und Minima
% Maxima = 1
% Minima = -1 
% Sattelpunkt = 0
psiI = zeros(1,length(phiI));
maxminI = zeros(1,length(phiI));
for i = 1:length(phiI)
    phi = phiI(i);
    H11 =  2 * (cos(2*phi)*(eyy-exx)-sin(2*phi)*gxy);
    H22 = -2 * (cos(phi)^2*exx+sin(phi)^2*eyy-ezz+0.5*sin(2*phi)*gxy);
    if H11*H22 > 0
        if H11 > 0
            maxminI(i) = -1;
        elseif H11 < 0
            maxminI(i) = 1;
        end
    end
end
EXTI = [phiI*180/pi;psiI*180/pi;maxminI];

%% Gleitung xz (ModeII)
phi = 0.5*atan(gxy/(exx-eyy));
phiII = [];
% finde alle Phi in [0 2pi]
while phi <= 2 * pi
    phiII = [phiII,phi];
    phi = phi + pi/2;
end
phiII = phiII(phiII > 0);
phiII = phiII(phiII <= 2* pi);
% Unterscheide Maxima und Minima
% Maxima = 1
% Minima = -1 
% Sattelpunkt = 0
psiII = repelem(pi/4,1,length(phiII));
maxminII = zeros(1,length(phiII));
for i = 1:length(phiII)
    phi = phiII(i);
    H11 =  2 * cos(2*phi)*(exx-eyy)+2*sin(2*phi)*gxy;
    H22 =  4 * (cos(phi)^2*exx+sin(phi)^2*eyy-ezz+0.5*sin(2*phi)*gxy);
    if H11*H22 > 0
        if H11 > 0
            maxminII(i) = -1;
        elseif H11 < 0
            maxminII(i) = 1;
        end
    end
end
EXTII = [phiII*180/pi;psiII*180/pi;maxminII];

%% Gleitung xy (ModeIII)
phi = 0.5*atan(-(exx-eyy)/gxy);
phiIII = [];
% finde alle Phi in [0 2pi]
while phi <= 2 * pi
    phiIII = [phiIII,phi];
    phi = phi + pi/2;
end
phiIII = phiIII(phiIII > 0);
phiIII = phiIII(phiIII <= 2* pi);
% Unterscheide Maxima und Minima
% Maxima = 1
% Minima = -1 
% Sattelpunkt = 0
psiIII = repelem(0,1,length(phiIII));
maxminIII = zeros(1,length(phiIII));
for i = 1:length(phiIII)
    phi = phiIII(i);
    H11 =  4 * sin(2*phi)*(exx-eyy)-4*cos(2*phi)*gxy;
    H22 =  sin(2*phi)*(exx-eyy)-cos(2*phi)*gxy;
    if H11*H22 > 0
        if H11 > 0
            maxminIII(i) = -1;
        elseif H11 < 0
            maxminIII(i) = 1;
        end
    end
end
EXTIII = [phiIII*180/pi;psiIII*180/pi;maxminIII];
end