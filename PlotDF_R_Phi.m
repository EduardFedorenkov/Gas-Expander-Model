function DF = PlotDF_R_Phi(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv, vz)
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp);
    VT0 = sqrt(2 * T0 / mg);
    Q = 1 / (sqrt(pi) * 2 * alpha^3) * diffCross * np ./ VTp.^2;
    Ueps = eps * VT0 / (2 * alpha) ./ VTp;
    nu = 16 * pi * alpha^3 * Q .* VTp.^3 .* exp(-Ueps.^2);
    % Make U grid, Ur and Uphi variables for plot 
    UGrid = linspace(-eps,eps,Nv);
    y0 = ones(1,Nv);
    Ur = y0'*UGrid;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;
    Uz = vz / VT0;
    f0 = n0 / (sqrt(pi) * VT0)^3 * exp( - UOrtSqr - Uz^2);
    DF = zeros(Nv, Nv, length(r));
    for i = 1:length(r)
        DF(:,:,i) = f0 .* exp(-nu(i) * a / VT0 * (r(i) / a * Ur ./ UOrtSqr  + sqrt(UOrtSqr - (r(i) / a)^2 * Uphi.^2) ./ UOrtSqr));
    end
    
    for i = 1:length(r)
        figure(i);
        pcolor(UGrid, UGrid, DF(:,:,i));
        shading flat;
        shading interp;
        title('F_{r\theta}(r = a)'); 
        xlabel('v_r/v_{Tg}'); 
        ylabel('v_\theta/v_{Tg}');
    end
end
