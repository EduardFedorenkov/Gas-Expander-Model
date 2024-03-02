function F = ComputeT0(n0, Trange, Twall, Pwall, kappa, np, Tp, mp, mg, diffCross, r, a, b, Nv)
    c = 3 * 10^10;
    eVtoErg = 1.6e-12;

    positiveIdxStart = fix(Nv/2) + 1;

    qr = zeros(1,length(Trange));
    Prr = zeros(1,length(Trange));
    for i = 1:length(Trange)
        [DFCold, gridStepCold, VrCold, VsqrCold] = PlotDFCold(np, Tp, mp, n0, Trange(i), mg, diffCross, r, a, Nv, false);
        nCold = Ncold(np, Tp, mp, n0, Trange(i), mg, diffCross, r, a, Nv);
        [DFhot, gridStepHot, VrHot, VsqrHot] = PlotDFHot(np, Tp, mp, Trange(i), nCold, mg, diffCross, r, a, Nv, false);
        
        PrrHot = 0;
        qrHot = 0;
        for j = 1:Nv
            PrrHot = PrrHot + sum(DFhot(:,positiveIdxStart:Nv,j,end) .* VrHot(:,positiveIdxStart:Nv).^2, "all");
            qrHot = qrHot + sum(DFhot(:,positiveIdxStart:Nv,j,end) .* VrHot(:,positiveIdxStart:Nv) .* VsqrHot(:,positiveIdxStart:Nv,j), "all");
        end
        PrrHot = PrrHot * mg * eVtoErg / c^2 * gridStepHot^3;
        qrHot = qrHot * mg * eVtoErg / c^2 * 0.5 * gridStepHot^3;
        
        PrrCold = 0;
        qrCold = 0;
        for j = 1:Nv
            PrrCold = PrrCold + sum(DFCold(:,positiveIdxStart:Nv,j,end) .* VrCold(:,positiveIdxStart:Nv).^2, "all");
            qrCold = qrCold + sum(DFCold(:,positiveIdxStart:Nv,j,end) .* VrCold(:,positiveIdxStart:Nv) .* VsqrCold(:,positiveIdxStart:Nv,j), "all");
        end
        PrrCold = PrrCold * mg * eVtoErg / c^2 * gridStepCold^3;
        qrCold = qrCold * mg * eVtoErg / c^2 * 0.5 * gridStepCold^3;
        qr(i) = qrHot + qrCold;
        Prr(i) = PrrHot + PrrCold;
    end
    F = (Trange * eVtoErg).^(3/2) + qr./Prr * Pwall / kappa * a *log(a/b) - (Twall * eVtoErg)^(3/2);
    plot(Trange,F);
end

