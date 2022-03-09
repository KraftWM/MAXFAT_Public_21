function plot_CRITPLANE(PHI,PSI,Z,dlflag,labelZ,varargin)
% Funktion erstellt 3D Plot aus Ergebnissen einer kritischen Ebenen
% Rechnung
%
% INPUT:
% PHI,PSI      - Winkel (Grad, nwinkel x 1)
%  Z           - Wert der auf die Z- Achse soll (nwinkel x 1)
% dlflag       - Z Werte sind Durchläufe   0  - nein 
%                                       sonst - Ja  
% labelZ       - Label Z-Achse
% varargin     - beliebige plot options
%__________________________________________________________________________

% ... Achensimits
maxPHI = max(PHI);
minPHI = min(PHI);
maxPSI = max(PSI);
minPSI = min(PSI);
if dlflag
    maxZ = 10^(ceil(log10(max(Z(Z<1e8)))));
    minZ = 10^(floor(log10(min(min(Z)))));
else
    maxZ = max(Z);
    minZ = min(Z);
end

% ... beschneide Z
Z(Z>1e8) = max(Z(Z<1e8));

% ... Anzahl Winkel
phi1 = PHI(1);
numphi = 1;
for i = 1 : length(PHI)
    if PHI(i) ~= phi1
        phi1 = PHI(i);
        numphi = numphi + 1;
    end
end

psi1 = PSI(1);
numpsi = 2;
while PSI(numpsi) ~= psi1
    numpsi = numpsi + 1;
end
numpsi = numpsi - 1;

% ... Erstelle Netzt
PHI = reshape(PHI,numpsi,numphi);
PSI = reshape(PSI,numpsi,numphi);
Z   = reshape(Z,numpsi,numphi);

% ... Erstelle Surface Plot
figure;
ax = gca;
latexinterpreter;
surf(ax,PHI,PSI,Z,varargin{:})
xlabel('$\varphi [{}^\circ]$');
ylabel('$\psi [{}^\circ]$');
zlabel(labelZ);

% ... Achensimits
if dlflag
    ax.ZScale = 'log';
end
ax.XLim = [minPHI, maxPHI];
ax.YLim = [minPSI, maxPSI];
ax.ZLim = [minZ, maxZ];


% ... Setzte View und colormap
shading(ax,'interp');
myjet = flipud(jet);
colormap(ax,myjet);
colorbar(ax);
ax.View = [0 90];
if dlflag
    set(ax,'ColorScale','log');
end


end % Ende Funktion