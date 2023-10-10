function intnCold = NcoldRtoL(nColdR, r, a, rk, sinPhiSqr, IsVrPositive, Nl)
    if(rk == 0)
        intnCold = trapz(r, nColdR);
    else
        cosPhi = sqrt(1 - sinPhiSqr);
        if(IsVrPositive == 1)
            lend = cosPhi + sqrt((rk/a)^2 - sinPhiSqr);
        else
            lend = cosPhi - sqrt((rk/a)^2 - sinPhiSqr);
        end
        l = linspace(0, lend, Nl);
        x = r/a;
        xint = sqrt(1 + l.^2 - 2*l*cosPhi);
        if(rk == a && xint(end) > 1)
            xint(end) = 1;
        end
        nColdL = interp1(x, nColdR, xint);
        intnCold = a * trapz(l, nColdL);
    end
end

