% run_AeroSim_Final.m
% A high-resolution, physics-based 3D aerodynamic analysis tool.
%
% --- FEATURES ---
% - Works directly with any 3D STL file.
% - Uses a robust hybrid physics model:
%   - Newtonian/Sine-Squared + Prandtl-Glauert for compression surfaces.
%   - Vacuum Pressure Limit (Cp_vac) for supersonic expansion.
%   - Compressibility-corrected skin friction for viscous drag.
% - Performs a highly detailed analysis across a wide flight envelope.
% - Fully implemented in the MKS unit system with results in kilonewtons (kN).
% - Automatically calculates reference area from the STL model.
% - All identified inaccuracies from previous versions have been corrected.

clear; clc; close all;

%% =================== 1. SIMULATION CONFIGURATION ===================
% --- Input File ---
stl_file ='put your stl file .stl'; % Place your STL file in the same directory

% --- High-Resolution Flight Envelope ---
% Defines a detailed matrix of points for a thorough analysis.
mach_subsonic   = 0.5:0.1:0.8;
mach_transonic  = 0.85:0.05:1.2; % Fine steps through the critical transonic region
mach_supersonic = 1.3:0.2:3.5;
mach_range = [mach_subsonic, mach_transonic, mach_supersonic];

altitude_range_m = 0:5000:25000; % Altitudes in METERS
aoa_range_deg = -2:1:12;         % Fine steps for angle of attack

%% =================== 2. SETUP AND INITIALIZATION ===================
fprintf('--- Starting High-Resolution 3D Aero-Analysis ---\n');

% Load STL geometry
fprintf('Loading STL file: %s\n', stl_file);
try
    model = stlread(stl_file);
catch ME
    error('Failed to read STL file. Ensure file exists and is a valid STL.');
end

% --- Automated Reference Area Calculation ---
S_ref_m2 = estimate_planform_area(model);
fprintf('Automated Planform Reference Area (S_ref): %.2f m^2\n', S_ref_m2);

% Pre-process geometry into panels
[panels] = preprocess_geometry(model);
fprintf('STL geometry processed into %d panels.\n', length(panels.areas));

% Initialize results storage
num_cases = length(mach_range) * length(altitude_range_m) * length(aoa_range_deg);
results = cell(num_cases, 1);
case_index = 1;

%% =================== 3. MAIN SIMULATION LOOP ===================
fprintf('Preparing to run simulation for %d detailed cases...\n', num_cases);
tic;

for alt_m = altitude_range_m
    % Get atmospheric conditions for this altitude (all in MKS units)
    [T, P, rho, mu] = get_standard_atmosphere(alt_m);
    a = sqrt(1.4 * 287 * T); % Speed of sound in m/s

    for M = mach_range
        % Freestream conditions in MKS units
        U_inf = M * a; % Freestream velocity in m/s
        q_inf = 0.5 * rho * U_inf^2; % Dynamic pressure in Pascals (N/m^2)
        
        L_char = sqrt(S_ref_m2); % Characteristic length for Reynolds number
        Re = rho * U_inf * L_char / mu;

        for aoa_deg = aoa_range_deg
            fprintf('  Running Case %d/%d: Alt=%.0f m, Mach=%.2f, AOA=%.1f deg\n', ...
                    case_index, num_cases, alt_m, M, aoa_deg);

            % --- Core Aerodynamic Calculation ---
            Cp = high_fidelity_solver(panels, M, aoa_deg);

            % --- Force Integration ---
            [L_N, D_N] = calculate_forces(panels, Cp, q_inf, aoa_deg, S_ref_m2, Re, M);

            % --- Store Results ---
            current_result.Altitude_m = alt_m;
            current_result.Mach = M;
            current_result.AOA_deg = aoa_deg;
            current_result.Lift_kN = L_N / 1000;   % Convert Newtons to kilonewtons
            current_result.Drag_kN = D_N / 1000;   % Convert Newtons to kilonewtons
            if D_N > 1e-6; current_result.LD_Ratio = L_N / D_N; else; current_result.LD_Ratio = inf; end
            
            results{case_index} = current_result;
            case_index = case_index + 1;
        end
    end
end
elapsed_time = toc;
fprintf('Simulation finished in %.2f seconds.\n', elapsed_time);

%% =================== 4. POST-PROCESSING AND VISUALIZATION ===================
fprintf('Generating results table and plots...\n');

% Convert results to a table for display
results_table = struct2table([results{:}]);

% --- Display Data Table ---
fig = figure('Name', 'Aerodynamic Performance Data (MKS Units)', 'NumberTitle', 'off', 'Position', [50 50 900 700]);
uitable(fig, 'Data', table2cell(results_table), 'ColumnName', results_table.Properties.VariableNames, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

% --- Generate Plots ---
visualize_results(results_table, model, panels, S_ref_m2);

fprintf('Analysis complete. Check the generated figures.\n');


%% =================== 5. GEOMETRY AND ATMOSPHERE FUNCTIONS ===================

function area = estimate_planform_area(model)
    % Estimates planform area by projecting vertices onto the XY plane.
    xy_points = model.Points(:, 1:2);
    try
        k = convhull(xy_points, 'Simplify', true);
        area = polyarea(xy_points(k,1), xy_points(k,2));
    catch
        % Fallback for co-linear or simple geometries
        area = (max(xy_points(:,1))-min(xy_points(:,1))) * (max(xy_points(:,2))-min(xy_points(:,2)));
    end
    if area < 1e-6; error('Failed to calculate a valid planform area. Check STL orientation.'); end
end

function [panel_data] = preprocess_geometry(model)
    % Processes the STL data into panels with centroids, areas, and normals.
    vertices = model.Points;
    faces = model.ConnectivityList;
    num_faces = size(faces, 1);
    
    panel_data.centroids = zeros(num_faces, 3);
    panel_data.normals = zeros(num_faces, 3);
    panel_data.areas = zeros(num_faces, 1);
    
    for i = 1:num_faces
        p1 = vertices(faces(i,1), :); p2 = vertices(faces(i,2), :); p3 = vertices(faces(i,3), :);
        panel_data.centroids(i,:) = (p1 + p2 + p3) / 3;
        cross_prod = cross(p2 - p1, p3 - p1);
        area = 0.5 * norm(cross_prod);
        panel_data.areas(i) = area;
        if area > 1e-9; panel_data.normals(i,:) = cross_prod / (2 * area); end
    end
end

function [T, P, rho, mu] = get_standard_atmosphere(h_m)
    % All inputs and outputs are in MKS units.
    g0 = 9.80665; R = 287.053; T0 = 288.15; P0 = 101325;
    layers = [0, T0, -0.0065; 11000, 216.65, 0; 20000, 216.65, 0.001; 
              32000, 228.65, 0.0028; 47000, 270.65, 0; 51000, 270.65, -0.0028;
              71000, 214.65, -0.002];
          
    idx = find(h_m >= layers(:,1), 1, 'last');
    h_b = layers(idx, 1); T_b = layers(idx, 2); L = layers(idx, 3);
    
    P_b = P0;
    if idx > 1
        for i = 2:idx
            h_prev = layers(i-1,1); T_prev=layers(i-1,2); L_prev=layers(i-1,3);
            if L_prev == 0; P_b = P_b * exp(-g0*(layers(i,1)-h_prev)/(R*T_prev));
            else; T_top = T_prev + L_prev*(layers(i,1)-h_prev); P_b = P_b * (T_top/T_prev)^(-g0/(L_prev*R)); end
        end
    end

    T = T_b + L * (h_m - h_b);
    if L == 0; P = P_b * exp(-g0 * (h_m - h_b) / (R * T));
    else; P = P_b * (T / T_b)^(-g0 / (L * R)); end
    rho = P / (R * T);
    
    mu_ref = 1.716e-5; T_ref = 273.15; S = 110.4;
    mu = mu_ref * (T / T_ref)^1.5 * (T_ref + S) / (T + S);
end

%% =================== 6. CORE PHYSICS AND SOLVER FUNCTIONS ===================

function Cp = high_fidelity_solver(panels, M, aoa_deg)
    % ** CORRECTED AND ROBUST SOLVER **
    gamma = 1.4;
    num_panels = length(panels.areas);
    normals = panels.normals;
    Cp = zeros(num_panels, 1);
    
    aoa_rad = deg2rad(aoa_deg);
    V_inf_dir = [cos(aoa_rad); 0; sin(aoa_rad)];
    
    dot_prod = normals * V_inf_dir;
    impact_indices = find(dot_prod < 0);
    shadow_indices = find(dot_prod >= 0);

    % --- 1. COMPRESSION SURFACES (Impact) ---
    delta = acos(-dot_prod(impact_indices));
    
    if M >= 1.0 % Supersonic Compression: Modified Newtonian Theory
        Cp(impact_indices) = 2 * sin(delta).^2;
    else % Subsonic Compression: Sine-Squared Law + Prandtl-Glauert
        Cp_inc = 2 * sin(delta).^2;
        beta = sqrt(1 - min(M^2, 0.9999)); % CORRECTED: Prevent M=1 singularity
        Cp(impact_indices) = Cp_inc / beta;
    end

    % --- 2. EXPANSION SURFACES (Shadow) ---
    if M >= 1.0 % Supersonic Expansion: Vacuum Pressure Limit
        % CORRECTED: More stable and physically robust than the previous P-M model.
        Cp_vac = -2 / (gamma * M^2);
        Cp(shadow_indices) = Cp_vac;
    else % Subsonic Expansion: Simple base pressure with P-G correction
        Cp(shadow_indices) = -0.1 / sqrt(1 - min(M^2, 0.9999));
    end
end

function [L_N, D_N] = calculate_forces(panels, Cp, q_inf, aoa_deg, S_ref, Re, M)
    % ** CORRECTED FORCE TRANSFORMATION AND FRICTION MODEL **
    
    % --- Pressure Forces ---
    force_magnitudes = Cp .* panels.areas * q_inf;
    force_vectors = -force_magnitudes .* panels.normals;
    F_total_pressure = sum(force_vectors, 1);
    
    % --- CORRECTED: Rotate forces from body-axis to wind-axis ---
    Fx_body = F_total_pressure(1);
    Fz_body = F_total_pressure(3);
    aoa_rad = deg2rad(aoa_deg);
    
    pressure_lift = Fz_body * cos(aoa_rad) - Fx_body * sin(aoa_rad);
    pressure_drag = Fz_body * sin(aoa_rad) + Fx_body * cos(aoa_rad);
    
    % --- Viscous (Skin Friction) Drag ---
    % CORRECTED: Use turbulent correlation with compressibility correction.
    Cf_inc = 0.074 / (Re^0.2); % Incompressible turbulent friction
    % Meador-Smart compressibility correction
    Cf = Cf_inc / (1 + 0.17 * M^2)^0.1295;
    
    wetted_area = sum(panels.areas);
    friction_drag = Cf * q_inf * wetted_area;
    
    % --- Total Lift and Drag in Newtons ---
    L_N = pressure_lift;
    D_N = pressure_drag + friction_drag;
    D_N = max(D_N, 0);
end

%% =================== 7. VISUALIZATION FUNCTION (CORRECTED) ===================

function visualize_results(results_table, model, panels, S_ref_m2)
    
    % --- 2D Plot: CL vs AOA for different Mach numbers ---
    figure('Name', 'Lift Coefficient vs. Angle of Attack', 'Position', [50, 500, 500, 400]);
    hold on; grid on;
    unique_mach = unique(results_table.Mach);
    colors = jet(length(unique_mach));
    
    % Plot for sea level (0 m)
    subset_sl = results_table(results_table.Altitude_m == 0, :);
    if isempty(subset_sl) && ~isempty(results_table)
        subset_sl = results_table(results_table.Altitude_m == min(results_table.Altitude_m),:);
        fprintf('Warning: No sea level data for CL plot. Using lowest altitude.\n');
    end

    if ~isempty(subset_sl)
        [~,~,rho_sl,~] = get_standard_atmosphere(subset_sl.Altitude_m(1));
        [T_sl,~,~,~] = get_standard_atmosphere(subset_sl.Altitude_m(1));
        a_sl = sqrt(1.4*287*T_sl);
        for i = 1:length(unique_mach)
            M = unique_mach(i);
            subset_mach = subset_sl(abs(subset_sl.Mach - M) < 1e-4, :);
            if ~isempty(subset_mach)
                q_inf = 0.5 * rho_sl * (M * a_sl)^2;
                CL = (subset_mach.Lift_kN * 1000) / (q_inf * S_ref_m2);
                plot(subset_mach.AOA_deg, CL, '-o', 'Color', colors(i,:), 'DisplayName', sprintf('M %.2f', M));
            end
        end
    end
    title('Lift Coefficient vs. AOA (at Lowest Altitude)');
    xlabel('Angle of Attack (deg)'); ylabel('Lift Coefficient (C_L)'); legend show;
    
    % --- 2D Plot: L/D Ratio vs. Mach Number ---
    figure('Name', 'L/D Ratio vs. Mach Number', 'Position', [550, 500, 500, 400]);
    hold on; grid on;
    unique_alt = unique(results_table.Altitude_m);
    colors = winter(length(unique_alt));
    % Plot for a constant AOA, e.g., 4 degrees
    aoa_plot = 4;
    if ~any(results_table.AOA_deg == aoa_plot) && ~isempty(results_table.AOA_deg)
        aoa_plot = results_table.AOA_deg(round(length(results_table.AOA_deg)/2));
        fprintf('Warning: AOA=4 not in range. Using AOA=%.1f for L/D plot.\n', aoa_plot);
    end
    subset_aoa = results_table(abs(results_table.AOA_deg - aoa_plot) < 1e-4, :);
    for i = 1:length(unique_alt)
       alt = unique_alt(i);
       subset_alt = subset_aoa(subset_aoa.Altitude_m == alt, :);
       if ~isempty(subset_alt)
           plot(subset_alt.Mach, subset_alt.LD_Ratio, '-s', 'Color', colors(i,:), 'DisplayName', sprintf('%.0f m', alt));
       end
    end
    title(sprintf('L/D Ratio vs. Mach Number (at AOA = %.0f deg)', aoa_plot));
    xlabel('Mach Number'); ylabel('Lift-to-Drag Ratio (L/D)'); legend show;

    % --- 3D Surface Plot: L/D Ratio vs. Mach and AOA ---
    figure('Name', 'L/D Performance Map', 'Position', [1050, 500, 500, 400]);
    if ~isempty(subset_sl) && height(subset_sl) > 3
        [X, Y] = meshgrid(unique(subset_sl.AOA_deg), unique(subset_sl.Mach));
        Z = griddata(subset_sl.AOA_deg, subset_sl.Mach, subset_sl.LD_Ratio, X, Y);
        if any(~isnan(Z(:)))
            surf(X, Y, Z, 'FaceColor','interp', 'EdgeAlpha', 0.2);
            title('L/D Ratio vs. Mach and AOA (at Lowest Altitude)');
            xlabel('Angle of Attack (deg)'); ylabel('Mach Number'); zlabel('L/D Ratio'); colorbar;
            view(30, 25);
        end
    end
    
    % --- Geometry Visualization: Pressure Distribution ---
    figure('Name', 'Pressure Distribution', 'Position', [50, 50, 800, 400]);
    
    % Select a high-supersonic case for visualization
    vis_case_idx = find(results_table.Mach >= 2.0 & abs(results_table.AOA_deg - 4) < 1e-4, 1);
    if isempty(vis_case_idx); vis_case_idx = find(abs(results_table.AOA_deg - 4) < 1e-4, 1, 'last'); end
    if isempty(vis_case_idx); vis_case_idx = 1; end % Failsafe
    
    vis_case = results_table(vis_case_idx,:);
    
    Cp_vis = high_fidelity_solver(panels, vis_case.Mach, vis_case.AOA_deg);
    
    % *** CORRECTED LINE HERE ***
    % Change 'interp' to 'flat'. This tells MATLAB to apply one color
    % per FACE, which matches the dimensions of our Cp_vis data.
    patch('Faces', model.ConnectivityList, 'Vertices', model.Points, ...
          'FaceVertexCData', Cp_vis, 'FaceColor', 'flat', 'EdgeColor', 'none');
          
    axis equal; view(3);
    camlight; lighting gouraud;
    title(sprintf('Pressure Coefficient (Cp) at M %.2f, AOA %.1f deg', vis_case.Mach, vis_case.AOA_deg));
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    colormap(jet); c = colorbar; c.Label.String = 'Pressure Coefficient (Cp)';
end
%% =================== 9. EXPORT RESULTS TO EXCEL SPREADSHEET ===================
% This section takes the final results table and writes it to a .xlsx file
% for easy analysis in Microsoft Excel or other spreadsheet software.

fprintf('Exporting full results table to Excel spreadsheet...\n');

% Define the output filename for the Excel file. It will be saved in the
% same directory as the script.
excel_filename = 'AeroSim_Final_Results.xlsx';

try
    % Use writetable, the modern and recommended function for this task.
    % It automatically writes the variable names from the table as headers.
    % 'WriteMode' 'replacefile' ensures that if the file exists, it's
    % overwritten with the new results.
    writetable(results_table, excel_filename, 'Sheet', 'Full Simulation Results', 'WriteMode', 'replacefile');
    
    % Confirm success to the user
    fprintf('Successfully exported %d data rows to the file "%s".\n', height(results_table), excel_filename);

catch ME
    % Catch potential errors. The most common error is that the user has
    % the Excel file open, which locks it and prevents MATLAB from writing to it.
    fprintf(2, '\n--- EXCEL EXPORT FAILED ---\n'); % fprintf(2,...) prints to the error stream (often in red)
    fprintf(2, 'An error occurred while trying to write to the Excel file.\n');
    fprintf(2, 'The most common reason for this is that the file "%s" is currently open in Excel.\n', excel_filename);
    fprintf(2, 'Please close the file and try running this section of the script again.\n\n');
    
    % Provide the specific MATLAB error for advanced debugging purposes.
    rethrow(ME);
end
%% =================== 8. L/D RATIO HEATMAP PLOTS BY ALTITUDE (CORRECTED) ===================
% This section generates a series of MATLAB figures, each displaying a
% heat map of the L/D ratio as a function of Mach number and Angle of
% Attack for a specific altitude.

fprintf('Generating L/D Ratio heat map plots for each altitude...\n');

% --- Determine Global Color Scale ---
% Using a global scale ensures that colors are comparable across all plots.
% A specific L/D value will have the same color regardless of altitude.
valid_ld_ratios = results_table.LD_Ratio(isfinite(results_table.LD_Ratio));
if ~isempty(valid_ld_ratios)
    global_ld_min = min(valid_ld_ratios);
    global_ld_max = max(valid_ld_ratios);
else
    global_ld_min = 0; % Fallback in case there are no valid L/D ratios
    global_ld_max = 1;
end

% --- Loop Through Each Altitude and Create a Plot ---
unique_altitudes = unique(results_table.Altitude_m);

for i = 1:length(unique_altitudes)
    current_alt = unique_altitudes(i);
    
    % Create a new figure for each altitude's heat map
    figure('Name', sprintf('L/D Heat Map @ %d m', current_alt), ...
           'Position', [100 + i*25, 100 + i*25, 650, 550], ...
           'NumberTitle', 'off');
    
    % Filter the main results table for data at the current altitude
    subset_alt = results_table(results_table.Altitude_m == current_alt, :);
    
    if height(subset_alt) < 3
        title(sprintf('Not enough data to generate heat map for %d m', current_alt));
        drawnow;
        continue;
    end
    
    % --- Pivot the data into a 2D grid for the heat map ---
    mach_grid = unique(subset_alt.Mach);
    aoa_grid = unique(subset_alt.AOA_deg);
    [X_grid, Y_grid] = meshgrid(aoa_grid, mach_grid);
    LD_grid = griddata(subset_alt.AOA_deg, subset_alt.Mach, subset_alt.LD_Ratio, X_grid, Y_grid);
    
    % --- Generate the Heat Map Plot ---
    h = heatmap(aoa_grid, mach_grid, LD_grid, 'Colormap', jet);
    
    % --- Customize the Plot Appearance ---
    h.Title = sprintf('L/D Ratio Performance Map at %d m', current_alt);
    h.XLabel = 'Angle of Attack (deg)';
    h.YLabel = 'Mach Number';
    h.ColorLimits = [global_ld_min, global_ld_max];
    h.MissingDataColor = [0.8, 0.8, 0.8];
    h.MissingDataLabel = 'No Data';
    h.CellLabelFormat = '%.2f';
    
    % --- CORRECTED SECTION ---
    % The HeatmapChart object 'h' does not have a .Colorbar property.
    % Instead, call the colorbar function to get a handle to the
    % colorbar object associated with the current axes.
    cb = colorbar;
    cb.Label.String = 'Lift-to-Drag Ratio (L/D)';
    % --- END CORRECTION ---
end

fprintf('Finished generating %d heat map figures.\n', length(unique_altitudes));
