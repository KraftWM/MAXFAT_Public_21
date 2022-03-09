function dTSED = neuberenergy(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS)
    dTSED = (SIG - REFSIG) .* dEPS + ... 
            (EPS - REFEPS) .* dSIG;
end