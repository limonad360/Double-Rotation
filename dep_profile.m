function Z = dep_profile(X, Y, center_x, center_y, C, val)

    %load('depline_Kaufman.mat');
    load('ExpData/depline_exp_130mm.mat');
    a = abs([center_x - max(X(1,:)), ...
                 center_x - min(X(1,:)), ...
                 center_y - max(Y(:,1)), ...
                 center_y - min(Y(:,1))]);
    
    r_max = max(a);
             
    r_max = r_max*1.1*sqrt(2);
    
    r = 0:(numel(depliney)+1)/profile_x_len:profile_x_len;
    
    T = depliney';

    ntt = 360;                             % // define how many angular division for the plot
    theta = linspace(0,2*pi,ntt);          % // create all the angular divisions
    [rr,tt] = meshgrid(r, theta);          % // generate a grid 

    z = repmat(T, ntt, 1);                 % // replicate "T" vector to match the grid
    [xr,yr,zr] = pol2cart(tt,rr,z);        % // convert everything to cartesian coordinates
    
    xr = xr + center_x;
    yr = yr + center_y;
    
    G = scatteredInterpolant(reshape(xr, [], 1), ...
                             reshape(yr, [], 1), ...
                             reshape(zr, [], 1), ...
                             'linear', 'nearest');
    
    Z = G(X, Y);
    
    if val == 1
        profile_norm = max(Z(:));
        Z = C * (Z./profile_norm);
    end
end