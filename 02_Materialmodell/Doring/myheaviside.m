function H = myheaviside(x)
% Implementierung der Heaviside Step Funktion weil Matlab version Irgendwie
% sehr langsam ist
H = sign(x);
H = 0.5 * (H + abs(H));
end
