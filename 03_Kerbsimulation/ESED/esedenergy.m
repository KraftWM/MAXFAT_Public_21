function dESED = esedenergy(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS)
    dESED = (SIG - REFSIG) .* dEPS;
end