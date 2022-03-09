%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Aitkens Delta Quadrat Verfahren                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dphat = deltaquadrat(dpiter,iter)
% Funktion zum beschleunigen der Konvergenz
    dphat = dpiter(iter) - ( dpiter(iter) - dpiter(iter-1) ).^2 / ...
        ( dpiter(iter) - 2*dpiter(iter-1) + dpiter(iter-2) );
    if dphat <= 0
        dphat = dpiter(iter);
    end
end