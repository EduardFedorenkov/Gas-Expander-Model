function F = ComputeT0DensePlasma(n0, Trange, Twall, Pwall, kappa, np, Tp, mp, mg, diffCross, a, b, Nv)
    c = 3 * 10^10;
    eVtoErg = 1.6e-12;
    kB = 1.38e-16;

    positiveIdxStart = fix(Nv/2) + 1;

    qr = zeros(1,length(Trange));
    Prr = zeros(1,length(Trange));
    for i = 1:length(Trange)
        beta23 = (2^(1/3) + 2^(-2/3));
        beta = beta23^(3/2);
        localCoef = 5^(2/3);
        nEff = 0.5 * n0 / (5 * beta) * (3/4) * ( sqrt(pi) * erf(localCoef * beta23) - ...
            2 * localCoef * beta23 * exp(-localCoef * beta23) );
        [DFhot, gridStepHot, VrHot, VsqrHot] = DFDensePlasma(np, Tp, mp, nEff, Trange(i), mg, diffCross, a, Nv, false);
        
        PrrHot = 0;
        qrHot = 0;
        for j = 1:Nv
            PrrHot = PrrHot + sum(DFhot(:,positiveIdxStart:Nv,j,end) .* VrHot(:,positiveIdxStart:Nv).^2, "all");
            qrHot = qrHot + sum(DFhot(:,positiveIdxStart:Nv,j,end) .* VrHot(:,positiveIdxStart:Nv) .* VsqrHot(:,positiveIdxStart:Nv,j), "all");
        end
        Prr(i) = PrrHot * mg * eVtoErg / c^2 * gridStepHot^3;
        qr(i) = qrHot * mg * eVtoErg / c^2 * 0.5 * gridStepHot^3;
        
    end
    F = (Trange * eVtoErg).^(3/2) + qr./Prr * Pwall / kappa * a *log(a/b) - (Twall * eVtoErg)^(3/2);
    plot(Trange,F);
end

