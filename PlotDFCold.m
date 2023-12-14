function [DF, gridStep, Vr, Vsqr] = PlotDFCold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv)
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
    gridStep = (UGrid(2) - UGrid(1)) * VT0;
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Vr = Ur * VT0;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;
    Vsqr = zeros(Nv, Nv, Nv);
    for i = 1:Nv
        Vsqr(:,:,i) = UOrtSqr(:,:) + UGrid(i)^2;
    end
    Vsqr = Vsqr * VT0^2;

    VelGridBoundary = eps * VT0;
    boundaryFactor = zeros(Nv,Nv,Nv);
    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) <= VelGridBoundary^2;
    end

    f0 = zeros(Nv, Nv, Nv);
    for i = 1:Nv
        f0(:,:,i) = n0 / (sqrt(pi) * VT0)^3 * exp( - UOrtSqr - Uz(i)^2);
    end
    DF = zeros(Nv, Nv, Nv, length(r));
    for j = 1:length(r)
        for i = 1:Nv
            DF(:,:,i,j) = f0(:,:,i) .* exp(-nu(j) * a / VT0 * (r(j) / a * Ur ./ UOrtSqr  + sqrt(UOrtSqr - (r(j) / a)^2 * Uphi.^2) ./ UOrtSqr));
        end
        DF(:,:,:,j) = DF(:,:,:,j) .* boundaryFactor;
    end
    
    for i = 1:length(r)
        figure(i);
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),i));
        shading flat;
        shading interp;
        title('F_{r\theta}(r = a)'); 
        xlabel('v_r/v_{Tg}'); 
        ylabel('v_\theta/v_{Tg}');
    end
end
