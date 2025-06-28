function [x, z] = panelgen(naca_code, n, aoa)

    %% ------------------------ enforce precond. ------------------------ %%
  
    arguments
        naca_code (1,4) char
        n (1,1) double = 100
        aoa (1,1) double = 0
    end
 
    assert(all(isstrprop(naca_code, 'digit')) && n > 0, 'naca_code must contain only digits and n must be positive.');

    %% ------------------------ parse naca_code ------------------------- %%

    % get diff. aerofoil params from naca_code
    m = str2double(naca_code(1)) * 1e-2;
    p = str2double(naca_code(2)) * 1e-1;
    t = str2double(naca_code(3:4)) * 1e-2;

    %% -------------------- calc. aerofoil geometry --------------------- %%
    x_coord = 1 - 0.5*(1 - cos(2*pi * (0:n)/n));

    % precomp. camber masks
    p_mask = x_coord < p;
    np_mask = ~p_mask;


    % def. seperate arr. for each idx to avoid repeated indexing
    x_less_p = x_coord(p_mask);
    x_greater_p = x_coord(np_mask);

    % calc. mean camber line
    y_camber = zeros(size(x_coord));
    y_camber(p_mask) = (m*p^-2) .* (2*p*x_less_p - x_less_p.^2);
    y_camber(np_mask) = (m*(1 - p)^-2) .* (1 - 2*p + 2*p*x_greater_p - x_greater_p.^2);

    % calc. thickness
    y_thickness = 5*t * (      ...
          0.2969*x_coord.^0.5  ...
        - 0.1260*x_coord       ...
        - 0.3516*x_coord.^2    ...
        + 0.2843*x_coord.^3    ...
        - 0.1015*x_coord.^4    ...
    );

    % calc. mean camber derivative
    y_camber_derivative = zeros(size(x_coord));
    y_camber_derivative(p_mask) = (p - x_less_p) .* (2*m*p^-2);
    y_camber_derivative(np_mask) = (p - x_greater_p) .* (2*m*(1 - p)^-2);

    % calc. mean camber angle
    theta = atan(y_camber_derivative);

    % get upper and lower ranges
    [~, leading_edge_idx] = min(x_coord);
    upper_idx = leading_edge_idx:length(x_coord);
    lower_idx = 1:leading_edge_idx;
    
    % init. coord arrays
    x = zeros(1, n+1);
    z = zeros(1, n+1);
    
    % calc. coords for upper and lower surfaces
    x(lower_idx) = x_coord(lower_idx) + y_thickness(lower_idx).*sin(theta(lower_idx));
    x(upper_idx) = x_coord(upper_idx) - y_thickness(upper_idx).*sin(theta(upper_idx));

    z(lower_idx) = y_camber(lower_idx) - y_thickness(lower_idx).*cos(theta(lower_idx));
    z(upper_idx) = y_camber(upper_idx) + y_thickness(upper_idx).*cos(theta(upper_idx));
    
    % merge trailing edge verts. to a single point
    x([1, end]) = 0.5*(x(1) + x(end));
    z([1, end]) = 0.5*(z(1) + z(end));
    
    % rotate wake panel
    wake_length = 1;
    wake_vec = complex(wake_length, 0) * exp(1i * -aoa); % why not ;)
    x = [x, x(end) + real(wake_vec)];
    z = [z, z(end)];

end
