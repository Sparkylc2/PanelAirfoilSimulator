%% -------------------------------- init. ------------------------------- %%
clear
clc

% the codebase is tiny so i find it's acceptable to use global vars. for
% ease of use

global NACA_CODE; % naca 4-series aerofoil code
global U_FS;      % freestream velocity
global AOA;       % angle of attack
global N;         % number of panels



global DOMAIN;            % plotting grid domain
global STREAM_X STREAM_Z; % streamline plotting grid (fine res.)
global QUIVER_X QUIVER_Z; % velocity field plotting grid (coarse res.)



% cache plotting domain and grids
DOMAIN = [-0.2, 1.2, -0.7, 0.7];
[STREAM_X, STREAM_Z] = meshgrid( ...
    linspace(DOMAIN(1), DOMAIN(2), 200), ...
    linspace(DOMAIN(3), DOMAIN(4), 200)  ...
);
[QUIVER_X, QUIVER_Z] = meshgrid( ...
    linspace(DOMAIN(1), DOMAIN(2), 20), ...
    linspace(DOMAIN(3), DOMAIN(4), 20)  ...
);


%% ----------------- get input & perform precond. chks. ----------------- %%


NACA_CODE = input('Enter your NACA-4 series aerofoil code (default: 2412): ', 's');
if isempty(NACA_CODE)
    NACA_CODE = '2412';
end

U_FS = input('Enter your freestream velocity (m/s) (default: 15): ');
if isempty(U_FS)
    U_FS = 15;
end
    
AOA = deg2rad(input('Enter your angle of attack (degrees) (default: 10): '));
if isempty(AOA)
    AOA = deg2rad(10);
end

N = input('Enter your number of panels (default: 200): ');
if isempty(N)
    N = 200;
end

 
%% --------------------------- xfoil import func. ----------------------- %%


function [aoa, cl] = xfoil_import(file_path)
    file = fopen(file_path, 'r');
    for i = 1:12 % skip header lines
        fgetl(file);
    end

    % get aoa and cl data
    data = textscan(file, '%f %f %*f %*f %*f %*f %*f', 'Delimiter', '', 'MultipleDelimsAsOne', true);
    fclose(file);

    aoa = deg2rad(data{1});
    cl  = data{2};
end


%% ------------------------- solve system func. ------------------------- %%

function mu = solve_system(panel_coord)
    % pull into scope
    global U_FS AOA N;

    % calc. preliminary vals. for panel geometry
    d = diff(panel_coord);
    beta = atan2(d(:, 2), d(:, 1));
    
    % precalc. req. panel coords.
    panel_midpoint = panel_coord(1:N, :) + 0.5*[d(1:N, 1), d(1:N, 2)]; 

    % precalc. trig funcs.
    cos_beta = cos(beta);
    sin_beta = sin(beta);

    % init. sys. matrices
    A = zeros(N + 1, N + 1); % lhs of sys.
    B = zeros(N + 1, 1); % rhs of sys.

    % apply kutta condition
    A(N + 1, [1, N, N + 1]) = [1, -1, 1];
    B(1:N) = -U_FS * sin(AOA - beta(1:N)); 
    
    % setup sys.
    for i = 1:N 
        for j = 1:N+1 % inc. last column (wake panel)
            [u_ij, v_ij] = cdoublet(panel_midpoint(i,:), panel_coord(j, :), panel_coord(j + 1, :)); % vel. at midpoint
            A(i, j) = v_ij*cos_beta(i) - u_ij*sin_beta(i); % perp. vel. component at midpoint
        end
    end
    mu = A\B; 
end

%% ---------------------------- velocity func. -------------------------- %%

function [u, v] = calculate_velocity(mesh_x, mesh_z, mu, panel_coord)
    % pull into scope
    global U_FS AOA N;
    
    % init. freestream vel. arr.
    dim = size(mesh_x);
    u = ones(dim) * (U_FS * cos(AOA));
    v = ones(dim) * (U_FS * sin(AOA));
    
    % set points inside to NaN
    outside = ~inpolygon(mesh_x, mesh_z, panel_coord(1:N, 1), panel_coord(1:N, 2));
    u(~outside) = NaN;
    v(~outside) = NaN;

    % get outside coords.
    points = [mesh_x(outside), mesh_z(outside)];
    
    % pre-alloc. arr. for induced vel.
    u_induced = zeros(dim);
    v_induced = zeros(dim);
    u_panel = zeros(N+1, 1);
    v_panel = zeros(N+1, 1);

    [r_idx, c_idx] = find(outside); % get outside idx.

    % iter. over points outside
    for i = 1:length(r_idx)
        % calc. induced vel. at each panel
        for k = 1:(N + 1)
            [u_panel(k), v_panel(k)] = cdoublet(points(i,:), panel_coord(k,:), panel_coord(k + 1, :));
        end
        
        % sum for point
        r_i = r_idx(i);
        c_i = c_idx(i);
        u_induced(r_i, c_i) = sum(mu .* u_panel);
        v_induced(r_i, c_i) = sum(mu .* v_panel);
    end

    % add induced vel. to freestream vel.
    u = u + u_induced;
    v = v + v_induced;
end

%% ---------------------------- main program ---------------------------- %%

% gen. aerofoil
[aerofoil_x, aerofoil_z] = panelgen(NACA_CODE, N, AOA);
aerofoil_coord = [aerofoil_x(:), aerofoil_z(:)];

% solve sys. and calc. vel. field
mu = solve_system(aerofoil_coord);
[stream_u, stream_v] = calculate_velocity(STREAM_X, STREAM_Z, mu, aerofoil_coord);
[quiver_u, quiver_v] = calculate_velocity(QUIVER_X, QUIVER_Z, mu, aerofoil_coord);

% calc. and print lift coeff.
cl = -2 * mu(N + 1) / U_FS;
fprintf("Lift Coefficient (C_L): %.3f\n", cl);

% plot aerofoil with streamlines
aerofoil_plot( ...
    aerofoil_coord, ...
    stream_u,       ...
    stream_v,       ...
    quiver_u,       ...
    quiver_v,       ...
    cl              ...
);

% ---------------------------- xfoil comparison -------------------------- %
if strcmp(NACA_CODE, '2412')

    % import xfoil data
    [aoa_xfoil, cl_xfoil] = xfoil_import('xf-naca2412-il-1000000.txt');

    % filter for range of interest
    filter = rad2deg(aoa_xfoil) >= 0  & rad2deg(aoa_xfoil) <= 10; 
    aoa_xfoil = aoa_xfoil(filter);
    cl_xfoil = cl_xfoil(filter);

    % init. iter. vars.
    N_Arr = [100, 200, 500];
    cl_foil = zeros(length(N_Arr), length(aoa_xfoil));
    mu = zeros(N + 1, 1); 


    U_FS = 15;  % update global state
    
    for i = 1:length(N_Arr)
        N = N_Arr(i); % update global state

        % gen. aerofoil 
        [aerofoil_x, aerofoil_z] = panelgen(NACA_CODE, N, 0);
        aerofoil_coord = [aerofoil_x(:), aerofoil_z(:)];

        % iter. over aoa vals. and solve sys. for each aoa to get cl
        for j = 1:length(aoa_xfoil)
            AOA = aoa_xfoil(j); % update global state

            % rotate wake panel
            wake_length = 1;
            wake_vec = complex(wake_length, 0) * exp(1i * -AOA);
            aerofoil_coord(end, :) = aerofoil_coord(end - 1, :) + [real(wake_vec), imag(wake_vec)];

            % solve. sys and calc. cl
            mu = solve_system(aerofoil_coord);
            cl_foil(i, j) = -2 * mu(N + 1) / U_FS;
        end
    end

    % plot last case
    [stream_u, stream_v] = calculate_velocity(STREAM_X, STREAM_Z, mu, aerofoil_coord);
    [quiver_u, quiver_v] = calculate_velocity(QUIVER_X, QUIVER_Z, mu, aerofoil_coord);
                
    aerofoil_plot( ...
        aerofoil_coord,  ... 
        stream_u,       ...
        stream_v,       ...
        quiver_u,       ...
        quiver_v,       ...
        cl_foil(i, end) ...
    );

    N = N_Arr; % update global state

    aoa_cl_plot( ...
        aoa_xfoil, ...
        cl_xfoil,  ...
        cl_foil    ...
    );
end



%% ------------------------ plotting funcs. ----------------------------- %%

% helper func. to ensure all plots have same label and grid config.
function label_plot(xlabel_text, ylabel_text, title_text)
    title(title_text, FontSize = 17, FontWeight = 'bold', Rotation = 0);
    xlabel(xlabel_text, FontSize = 15, FontWeight = 'light', Rotation = 0);
    ylabel(ylabel_text, FontSize = 15, FontWeight = 'light', Rotation = 0);
    grid on;
end

% helper callback func. to reposition annotation dynamically rel. to ax. 
function adjust_annotation(ann, ax)
    x_limits = xlim(ax);
    y_limits = ylim(ax);

    % recalc. annotation position based on new axis limits
    new_x = x_limits(2) + 0.05 * (x_limits(2) - x_limits(1)); 
    new_y = y_limits(1) + 0.75 * (y_limits(2) - y_limits(1));

    % update annotation position
    set(ann, Position = [new_x, new_y]);
end

% helper func. to add annotation to aerofoil plots
function create_annotation(fig, ax, annotation_text)
    % create annotation box
    ann = text( ...
        0,                           ...
        0,                           ...
        annotation_text,             ...
        FontSize = 14,               ...
        FontName = 'Helvetica Neue', ...
        FontWeight = 'light',        ...
        Color = 'k',                 ...
        EdgeColor = 'k',             ...
        BackgroundColor = 'w',       ...
        LineWidth = 1,               ...
        Margin = 5,                  ...
        VerticalAlignment = 'top',   ...
        HorizontalAlignment = 'left' ...
    );
    % set initial pos.
    adjust_annotation(ann, ax);
    % set callback func.
    set(fig, 'SizeChangedFcn', @(src, event) adjust_annotation(ann, ax));
end



function aerofoil_plot(aerofoil_coord, stream_u, stream_v, quiver_u, quiver_v, cl)
    % pull into scope
    global NACA_CODE N AOA U_FS DOMAIN STREAM_X STREAM_Z QUIVER_X QUIVER_Z;
    

    % update aerofoil coord. to exclude wake panel for plotting
    aerofoil_coord = aerofoil_coord(1:N, :);
    % ------------------------ streamline plot ------------------------ %
    % gen. fig., grab ax. ref., and shift left to give space for annotation
    fig1 = figure(Name = sprintf('Streamlines - NACA %s', NACA_CODE));
    ax1 = gca;
    set(ax1, Position = [0.08, 0.15, 0.65, 0.75]); 

    % plot aerofoil (without wake panel)
    plot( ...
        aerofoil_coord(:, 1), ...
        aerofoil_coord(:, 2), ...
        Color = 'k',          ...
        LineWidth = 2         ...
    );

    hold on;

    % plot streamlines
    s = streamslice( ...
        STREAM_X, ...
        STREAM_Z, ...
        stream_u, ...
        stream_v  ...
    );
    set(s, Color = 'r', LineWidth = 1);
    
    % set plot props.
    label_plot('x/c', 'z/c', sprintf('Streamlines around NACA %s', NACA_CODE));
    axis(DOMAIN);
    pbaspect([(DOMAIN(2)-DOMAIN(1)) (DOMAIN(4)-DOMAIN(3)) 1]);

    hold off;

    % add annotation
    annotation_text = sprintf( ...
        '\\alpha = %.1f°\nU_\\infty = %.1f (m/s)\nN = %d\nC_L = %.3f', ...
        rad2deg(AOA), U_FS, N, cl                                      ...
    );
    create_annotation(fig1, ax1, annotation_text);

    saveas(fig1, sprintf('NACA_%s_streamlines_AoA_%.1f_N_%d_Ufs_%.1f.png', NACA_CODE, rad2deg(AOA), N, U_FS));
    

    % ---------------------- velocity field plot ---------------------- %
    % gen. fig., grab ax. ref., and shift left to give space for annotation
    fig2 = figure(Name = sprintf('Velocity Vector Field - NACA %s', NACA_CODE));
    ax2 = gca;
    set(ax2, Position = [0.08, 0.15, 0.65, 0.75]); 

    % plot aerofoil (without wake panel)
    plot( ...
        aerofoil_coord(:, 1), ...
        aerofoil_coord(:, 2), ...
        Color = 'k',          ...
        LineWidth = 2         ...
    );

    hold on;

    % plot velocity vec. field
    quiver( ...
        QUIVER_X,      ...
        QUIVER_Z,      ...
        quiver_u,      ...
        quiver_v,      ...
        LineWidth = 1, ...
        Color = 'b'    ...
    );


    % set plot props.
    label_plot('x/c', 'z/c', sprintf('Velocity field around NACA %s', NACA_CODE));
    axis(DOMAIN);
    pbaspect([(DOMAIN(2)-DOMAIN(1)) (DOMAIN(4)-DOMAIN(3)) 1]);

    hold off;

    % add annotation
    annotation_text = sprintf( ...
        '\\alpha = %.1f°\nU_\\infty = %.1f (m/s)\nN = %d\nC_L = %.3f', ...
        rad2deg(AOA), U_FS, N, cl                                      ...
    );
    create_annotation(fig2, ax2, annotation_text);


    saveas(fig2, sprintf('NACA_%s_arrows_AoA_%.1f_N_%d_Ufs_%.1f.png', NACA_CODE, rad2deg(AOA), N, U_FS));
end


function aoa_cl_plot(aoa_xfoil, cl_xfoil, cl_foil)
    % pull into scope
    global NACA_CODE U_FS N
    
    % grab fig. ref.
    fig = figure(Name = sprintf('NACA %s Lift Curve', NACA_CODE));

    % init. plot props
    markers = ["square", "diamond", "^"];
    colors = ['y', 'g', 'b'];

    % plot cl vs aoa for xfoil
    plot( ...
        rad2deg(aoa_xfoil),    ...
        cl_xfoil,              ...
        LineStyle = '-',       ...
        LineWidth = 1.5,       ...
        Color = 'k',           ...
        MarkerSize = 6,        ...
        Marker = 'o',          ...
        MarkerFaceColor = 'k', ...
        MarkerEdgeColor = 'k', ...
        DisplayName = 'XFOIL'  ...
    );

    hold on;

    % iter. over N values and plot cl vs aoa for panel method
    for i = 1:length(N)
        plot( ...
            rad2deg(aoa_xfoil),                    ...
            cl_foil(i, :),                         ...
            LineStyle = '-',                       ...
            LineWidth = 1.5,                       ...
            Color = colors(i),                     ...
            MarkerSize = 6,                        ...
            Marker = markers(i),                   ...
            MarkerFaceColor = colors(i),           ...
            MarkerEdgeColor = 'k',                 ...
            DisplayName = sprintf('N = %d', N(i))  ...
        );
    end

    hold off

    % set axis props.
    axis auto;
    label_plot("α", "C_L", sprintf('NACA %s Lift Curve at U_∞ = %.1f (m/s)', NACA_CODE, U_FS));   
    legend('show', Location = 'northwest');

    saveas(fig, sprintf('NACA_%s_liftcurve_Ufs_%.1f.png', NACA_CODE, U_FS));
end
