function dPHI = modneuberenergy(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS,...
                                EEPS,REFEEPS,dEEPS)
% obacht rückgabewert keine Energie
    dTSED = (SIG - REFSIG) .* dEPS + ... 
            (EPS - REFEPS) .* dSIG;
%     phi2 = (EEPS(2) - REFEEPS(2)) * dEPS(1) - ...        % ee_zz*de_yy
%            (EEPS(1) - REFEEPS(1)) * dEPS(2) + ...        % ee_yy*de_zz
%            ( EPS(1) -  REFEPS(1)) * dEEPS(2) - ...       % e_yy*dee_zz 
%            ( EPS(2) -  REFEPS(2)) * dEEPS(1);            % e_zz*dee_yy 
%     phi2 = (EEPS(2) - REFEEPS(2)) * dEPS(1) - ...        % ee_zz*de_yy
%            (EEPS(1) - REFEEPS(1)) * dEPS(2);             % ee_yy*de_zz
    phi2 = ( EPS(1) -  REFEPS(1)) * dEEPS(2) - ...       % e_yy*dee_zz 
           ( EPS(2) -  REFEPS(2)) * dEEPS(1);            % e_zz*dee_yy            
    dPHI = [dTSED(1) + dTSED(2); ...
                           phi2; ...
                      dTSED(3)];
end