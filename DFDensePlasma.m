function [DF, gridStep, Vr, Vsqr] = DFDensePlasma(np, Tp, mp, nEff, T0, mg, diffCross, a, Nv, isPlot)
    c = 3 * 10^10;
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp) * c;
    VT0 = sqrt(2 * T0 / mg) * c;
    Q = 1 / (sqrt(pi) * 2 * alpha^3) * diffCross * np / VTp^2;
    Ueps = eps * VT0 / (2 * alpha) / VTp;
    nu = 16 * pi * alpha^3 * Q * VTp^3 * exp(-Ueps^2);
    da = 5 * VT0 / nu;

    if (a * 0.3 < da)
        disp('Exeption: da seams too large!!!');
    end

    P = zeros(Nv, Nv, Nv);
    DF = zeros(Nv, Nv, Nv);

    VelGridBoundary = eps * VT0 / VTp;                       % V < Veps domen is no need for this solution
    boundaryFactor = zeros(Nv,Nv,Nv);

    UGrid = linspace(-eps,eps,Nv);
    gridStep = (UGrid(2) - UGrid(1)) * VTp;
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Vr = Ur * VTp;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;
    UOrtSqrt = sqrt(UOrtSqr);
    Vsqr = zeros(Nv, Nv, Nv);
    for i = 1:Nv
        Vsqr(:,:,i) = UOrtSqr(:,:) + UGrid(i)^2;
    end
    Vsqr = Vsqr * VTp^2;

    lenTypeOne = 2 * a * abs(Ur) ./ UOrtSqrt .* ((abs(Uphi) ./ UOrtSqrt) >= (1 - da / a));
    lenTypeTwo = 2 * a * abs(Ur) ./ UOrtSqrt .* (1 - sqrt(1 - (UOrtSqr ./ Ur.^2) * (da/a)^2 * (2*a/da -1))) ...
        .* ((abs(Uphi) ./ UOrtSqrt) < (1 - da / a));

    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) >= VelGridBoundary^2;
        P(:,:,i) = (2 * alpha) * Q ./ sqrt(UOrtSqr + Uz(i)^2) .* exp(- (UOrtSqr + Uz(i)^2) ./ (4 * alpha^2));
    end
    P = P .* boundaryFactor;

    for i = 1:Nv
        DF(:,fix(Nv/2):Nv,i) = P(:,fix(Nv/2):Nv,i) ./ UOrtSqrt(:,fix(Nv/2):Nv) .* ...
            ( lenTypeOne(:,fix(Nv/2):Nv) + lenTypeTwo(:,fix(Nv/2):Nv) );
    end
    DF = DF * (nEff / VTp);

    if (isPlot)
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2)));
        shading flat;
        shading interp;
        title('f_g(r = 0)'); 
        xlabel('v_r/v_{Tp}'); 
        ylabel('v_\phi/v_{Tp}');
        colormap('jet');
        colorbar;
        set(gca,'ColorScale','log');
        hold on
        contour(UGrid, UGrid, DF(:,:,fix(Nv/2)), [1e-12 1e-9]);
        hold off
    end
end

