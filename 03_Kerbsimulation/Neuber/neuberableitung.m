function [G, P] = neuberableitung(SIG,EPS,...
                                  REFSIG,REFEPS)
% =====================================================================
% P
P = diag(SIG - REFSIG);

% =====================================================================
% G
G = diag(EPS - REFEPS);

end % Ende Funktion