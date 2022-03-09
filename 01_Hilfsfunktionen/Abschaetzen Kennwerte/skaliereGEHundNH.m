function q = skaliereGEHundNH(fwt)
% q(x) = a - b * f(x)
% q(1/sqrt(3)) = 0    -> GEH
% q(1)         = 1    -> NH
% f(1/sqrt(3)) = f13
% f(1)         = f1
% a =  f13/(f13-f1)
% b = 1/(f13-f1)
% Hier Ansätze  f(x) = (1/x)^d (d=1 entspricht RiLiNiLi)
% -------------------------------------------------------------
d = 1;
f1 = 1;
f13 = sqrt(3)^d;
q = (f13 - 1./fwt.^d)./(f13-f1);
end