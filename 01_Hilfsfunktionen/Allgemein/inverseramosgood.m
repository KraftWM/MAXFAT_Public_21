function SIG = inverseramosgood(EPS,E,K,n)
% Invertiere Ramberg Osgood Gleichung
% eps = sig/E + (sig/K)^(1/n)

SIG = zeros(size(EPS));

for i = 1:length(EPS)
    eps = EPS(i);
    % ... Startwert
    if eps > 0.0005
        sig = K*eps^n;
    else
        sig = eps*E;
    end

    f = zielfun(eps,sig,E,K,n);
    tol = 1e-14;
    iter = 1;
    dsig = 1;
    while abs(f) > tol
        if iter > 10000
            error('keine Konvergenz')
        end
        f1 = zielfun(eps,sig+dsig,E,K,n);
        f2 = zielfun(eps,sig-dsig,E,K,n);
        df = (f1-f2)/(2*dsig);
        sig = sig - f/df;
        f = zielfun(eps,sig,E,K,n);
        iter = iter + 1;
    end
    SIG(i) = sig;
end

end % Ende Funktion

function f = zielfun(eps,sig,E,K,n)
    f = eps - sig./E - (sig/K)^(1/n);
end
