function nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv)
    c = 3 * 10^10;
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp) * c;
    VT0 = sqrt(2 * T0 / mg) * c;
    Q = 1 / (sqrt(pi) * 2 * alpha^3) * diffCross * np ./ VTp.^2;
    Ueps = eps * VT0 / (2 * alpha) ./ VTp;
    nu = 16 * pi * alpha^3 * Q .* VTp.^3 .* exp(-Ueps.^2);
    % Make U grid, Ur and Uphi variables for plot 
    UGrid = linspace(-eps,eps,Nv);
    gridStep = UGrid(2) - UGrid(1);
    y0 = ones(1,Nv);
    Ur = y0'*UGrid;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;
    Uz = UGrid;
    f0 = zeros(Nv, Nv, Nv);
    factor = zeros(Nv, Nv, Nv);
    DF = zeros(Nv, Nv, Nv, length(r));
    nCold = zeros(size(r));
    for z = 1 : Nv
        factor(:,:,z) = (UOrtSqr + Uz(z)^2) <= eps^2;
        f0(:,:,z) = n0 / (sqrt(pi) * VT0)^3 * exp( - UOrtSqr - Uz(z)^2);
    end
    
    for i = 1:length(r)
        for z = 1 : Nv 
            DF(:,:,z,i) = factor(:,:,z) .* f0(:,:,z) .* exp(-nu(i) * a / VT0 * (r(i) / a * Ur ./ UOrtSqr  + sqrt(UOrtSqr - (r(i) / a)^2 * Uphi.^2) ./ UOrtSqr));
        end
    end

    for i = 1:length(r)
        nCold(i) = sum(DF(:,:,:,i), "all");
    end
    nCold = nCold * (gridStep * VT0)^3;
end

