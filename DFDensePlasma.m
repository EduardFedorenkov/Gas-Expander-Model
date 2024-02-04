function [DF, gridStep, Vr, Vsqr] = DFDensePlasma(np, Tp, mp, nEff, T0, mg, diffCross, Nv)
    c = 3 * 10^10;
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp) * c;
    VT0 = sqrt(2 * T0 / mg) * c;
    Q = 1 / (sqrt(pi) * 2 * alpha^3) * diffCross * np / VTp^2;
    Ueps = eps * VT0 / (2 * alpha) / VTp;
    nu = 16 * pi * alpha^3 * Q * VTp^3 * exp(-Ueps^2);
    da = 5 * VT0 / nu;

    P = zeros(Nv, Nv, Nv);
    DF = zeros(Nv, Nv, Nv);

    VelGridBoundary = eps * VT0 / max(VTp);                       % V < Veps domen is no need for this solution
    boundaryFactor = zeros(Nv,Nv,Nv);

    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) >= VelGridBoundary^2;
        P(:,:,i) = (2 * alpha) * Q / sqrt(UOrtSqr + Uz(i)^2) .* exp(- (UOrtSqr + Uz(i)^2) ./ (4 * alpha^2));
    end
    P = P .* boundaryFactor;

    
end

