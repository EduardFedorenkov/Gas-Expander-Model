function DF = PlotDFHot(np, Tp, mp, T0, Ncold, mg, diffCross, r, a, Nv)
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp);
    VT0 = sqrt(2 * T0 / mg);
    Q = diffCross ./ (2 * alpha^3 * sqrt(pi) * VTp.^2) .* np;
    VelGridBoundary = eps * VT0 / max(VTp);                       % V < Veps domen is no need for this solution

    UGrid = linspace(-eps,eps,Nv);
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;

    P = zeros(Nv, Nv, Nv,length(r));
    DF0 = zeros(Nv, Nv, Nv);
    DF = zeros(Nv, Nv, Nv, length(r)-1);
    boundaryFactor = zeros(Nv,Nv,Nv);
    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) >= VelGridBoundary;
    end

    intNcold0 = trapz(r,Ncold);
    intNcold = zeros(Nv, Nv, length(r));
    x = (r / a);
    for k = 2:length(r)
        for j = 1:Nv
            for i = 1:Nv
                sinPhiSqr = (r(k) / a)^2 * UGrid(i)^2 / (UGrid(j)^2 + UGrid(i)^2);
                if UGrid(j) <= 0
                    integrationDomen = (x.^2 > max([sinPhiSqr, x(k)^2])) ./ sqrt(1 - sinPhiSqr ./ x.^2);
                    intNcold(i, j, k) = trapz(x, integrationDomen .* Ncold);
                else
                    integrationDomen = (x.^2 > sinPhiSqr) ./ sqrt(1 - sinPhiSqr ./ x.^2);
                    intNcold(i, j, k) = trapz(x, integrationDomen .* Ncold);
                    if x(k)^2 > sinPhiSqr
                        integrationDomen2 = (x <= x(k)) .* integrationDomen;
                        intNcold(i, j, k) = intNcold(i, j, k) + ...
                            trapz(x, integrationDomen2 .* Ncold);
                    end
                end
            end
        end
    end

    for j = 1:length(r)
        for i = 1:Nv
            P(:,:,i,j) = (2 * alpha) * Q(j) ./ sqrt(UOrtSqr + Uz(i)^2) .* exp(- (UOrtSqr + Uz(i)^2) ./ (4 * alpha^2));
        end
        P(:,:,:,j) = P(:,:,:,j) .* boundaryFactor;
    end

    for i = 1:Nv
        DF0(:,:,i) = P(:,:,i,1) ./ sqrt(UOrtSqr) * intNcold0;
    end
    DF0 = DF0 * a;

    DF(:,:,:,1) = DF0;
    for k = 2:length(r)
        for j = 1:Nv
            DF(:,:,j,k) = P(:,:,j,k) ./ sqrt(UOrtSqr(:,:)) .* intNcold(:,:,k);
        end
    end
    DF = DF * a;
    
    figure(1);
    pcolor(UGrid, UGrid, DF0(:,:,fix(Nv/2)));
    shading flat;
    shading interp;
    title('F(r = 0)'); 
    xlabel('v_r/v_{Tp}'); 
    ylabel('v_\phi/v_{Tp}');

    for i = 2:length(r)
        figure(i);
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),i));
        shading flat;
        shading interp;
        title('F(r = a)'); 
        xlabel('v_r/v_{Tp}'); 
        ylabel('v_\phi/v_{Tp}');
    end

end

