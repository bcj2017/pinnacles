% Initialize spatial fields
function [u,v,T,phi] = initialize(rr,zz,Nr,Nz,Tstar,T_inf,H,L,pinnacleSize,perturbationVal,perturbationHeightCm)
    % New nondimensionalized temperatures:
    theta_inf = (T_inf - Tstar)/T_inf;
    theta_0 = -Tstar / T_inf;

    % Set up flow, temp, phase fields
    u = zeros(Nr-1,Nz);
    v = zeros(Nr,Nz-1);
    T = theta_inf*ones(Nr,Nz);
    phi = ones(Nr,Nz);
    
    % Standard equilibrium shape
    R = 1/L * pinnacleSize * 1e-6; % middle number is microns %R = .8/.17 * 6e-4; %R = 1e-2;
    t0 = asin(sqrt((R-sqrt(R^2 - 4*(R/4-H)*(3*R/4)))/(2*(R/4-H)))); %final tangent angle
    % if perturbationVal ~= 0
    theta = linspace(pi/2,t0,1e3);
    % else
        % theta = linspace(pi/2,t0);
    % end
    Xr = R*(cos(theta)./sin(theta).^3); %r coordinates
    Xz = H-R*(-1./sin(theta).^2 + (3/4)./sin(theta).^4 + (1/4)); %z coordinates

    % NEW: reinterpolate equispaced in x (useful for interpolating new
    % points)
    nnpt = 1e3; length(theta);1e3; % theta spacing gives bad interpolation for 100 pts!
    Xrinterp = linspace(min(Xr),max(Xr),nnpt);
    Xz = interp1(Xr,Xz,Xrinterp);
    Xr = Xrinterp;
    dx = Xr(2) - Xr(1); % x spacing for finite difference
        
    % Perturbation of cm height perturbationHeightCm
    pertLengthDimless = perturbationHeightCm * 1e-2 * 1/L; % convert to meters and nondimensionalize
    if perturbationVal == 1 % sharp
        Hmax = max(Xz) + pertLengthDimless; % Identify new max height
        % Identify points outside of 0.5cm of original top
        indices = find(Xz < (max(Xz) - pertLengthDimless)); % to keep
        apexIndices = find(Xz >= (max(Xz) - pertLengthDimless)); % to delete
    elseif perturbationVal == -1 % blunt
        Hmax = max(Xz) - pertLengthDimless;
        % Identify points within 0.5cm below *new* top
        indices = find(Xz < (Hmax - pertLengthDimless)); % to keep
        apexIndices = find(Xz >= (Hmax - pertLengthDimless)); % to delete
    end
    % Preserve original pinnacle
    Xr_unperturbed = Xr; Xz_unperturbed = Xz; 
    % Perturbations given an Hmax and points to delete
    if perturbationVal ~= 0 % yes to perturbations
        % Delete points that we will overwrite
        Xr = Xr_unperturbed(indices); Xz = Xz_unperturbed(indices);
        Xr_deletedApex = Xr_unperturbed(apexIndices); Xz_deletedApex = Xz_unperturbed(apexIndices);
        
        % Interpolation points
        xApex = 0; zApex = Hmax; % pre-determined new apex point
        xBody = min(Xr); zBody = max(Xz); % end point of cropped body

        % Interpolation slope for apex: zero slope at apex
        mApex = 0; % slope of zero at the apex itself
        % Interpolation slope along body: compute mBody using a centered finite
        % difference f'(x) = (f(x+h) - f(x-h)) / (2h)
        xMinusH = max(Xr_deletedApex); zMinusH = min(Xz_deletedApex); % left point
        % Find second-smallest point in Xr and its index
        [vals, idx] = unique(Xr, 'sorted'); secondVal = vals(2); secondIdx = find(Xr == secondVal, 1);
        xPlusH = secondVal; zPlusH = Xz(secondIdx);
        % Use these points to get slope along body:
        mBody = (zPlusH-zMinusH)/(2*dx);


        % Cubic Hermite spline (equation from Wikipedia)
        % Step 1: interpolation x values
        xInterp = Xr_deletedApex;
        % Confirm that xInterp includes zero
        disp(['Minimum x value in deleted apex: ',num2str(min(Xr_deletedApex))])
        % xInterp(1) = 0; % 1e-20 but whatever, not sure if this is useful DELETE PROBS
        % Step 2: set variable t for change of variables
        t = (xInterp - xApex)/(xBody - xApex);
        % Step 3: define Hermite basis functions
        h00 = @(tt) 2*tt.^3 - 3*tt.^2 + 1;
        h10 = @(tt) tt.^3 - 2*tt.^2 + tt;
        h01 = @(tt) -2*tt.^3 + 3*tt.^2;
        h11 = @(tt) tt.^3 - tt.^2;
        % Step 4: Calculate interpolated z values using the Hermite spline
        zInterp = h00(t) * zApex + h10(t) * mApex * (xBody - xApex) + ...
                  h01(t) * zBody + h11(t) * mBody * (xBody - xApex);
        % Finished with cubic spline!

        % Frankenstein together interpolation and spline
        Xr = [xInterp Xr];
        Xz = [zInterp Xz];
    end

    % Extend Xr and Xz, determine ice boundary, set temp and phase
    % accordingly
    Xr = [0,0,Xr]; Xz = [0,0,Xz];
    inds = inpolygon(rr,zz,Xr,Xz);
    T(inds) = theta_0; phi(inds) = 0; % temperature and phase in solid


    
end
