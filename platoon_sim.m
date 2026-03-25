clear;
clc;
close all;

%% ================== 1. User Configuration ==================
% --- Select simulation scenario (1: Low penetration Type 1, 2: High penetration Type 2) ---
SCENARIO_ID = 1; 

% --- Key time parameters ---
deta_t = 0.1;          % time step [s]
t_init = 2;            % initialization time (all vehicles at equilibrium) [s]
t_stable_run = 48;     % time after initialization before perturbation (leader cruises) [s]
t_post_stable = 50;    % time after perturbation to stabilize [s]

% --- Basic simulation settings ---
Num_platoon = 10;       % total number of vehicles (excluding leader)
Ve = 20;                % equilibrium speed [m/s]

%% ================== 2. Load Data and General Parameters ==================
try
    load acc_ngsim.mat
    % Read perturbation data (leader deceleration profile)
    HDV_deacceleration = acceleration(4,1:326);
catch
    warning('Data file not found. Using default simulated data.');
    HDV_deacceleration = -2 * ones(1, 326); 
end

% Total simulation time and steps
% Phase 1: initialization (t_init)
% Phase 2: steady run (t_stable_run)
% Phase 3: perturbation (length of HDV_deacceleration)
% Phase 4: post-perturbation stabilization (t_post_stable)
T = t_init + t_stable_run + length(HDV_deacceleration)*deta_t + t_post_stable;
total_steps = round(T/deta_t);

% --- Vehicle physical limits ---
max_acc = 5;    % [m/s²]
min_acc = -5;   % [m/s²]
S_go = 50;      % [m] gap at which desired speed becomes V_max
S_st = 5;       % [m] standstill gap
V_max = 30;     % [m/s] maximum speed

% --- Driver/vehicle parameters (common) ---
h     = 1.2;            % CAV desired time gap [s]
theta = round(0.1/deta_t); % communication delay [steps]
tao   = 0.5;            % vehicle time lag [s]
D     = 5;              % standstill spacing [m]

% --- HDV heterogeneity parameters (OVM model) ---
alpha = repmat(0.2, 1, Num_platoon);  
beta  = repmat(0.4, 1, Num_platoon);
th    = repmat(1.5, 1, Num_platoon);  
phi   = repmat(10, 1, Num_platoon);   

%% ================== 3. Scenario Configuration ==================
% --- 3.1 Vehicle platoon composition based on scenario ---
switch SCENARIO_ID
    case 1 % Type 1 (low penetration)
        Vehicle_type = [0 0 1 0 0 1 0 0 0 1];   % 0: HDV, 1: CAV
        Num_HDV      = [0 1 2 0 1 2 0 1 2 3];    % number of HDVs ahead of each CAV
        index        = [0 0 1 0 0 2 0 0 0 3];    % index for CACCm gains (position in platoon)
        
    case 2 % Type 2 (high penetration)
        Vehicle_type = [0 1 1 0 1 1 0 1 1 1]; 
        Num_HDV      = [0 1 0 0 1 0 0 1 0 0]; 
        index        = [0 1 0 0 2 0 0 3 0 0]; 
        
    otherwise
        error('Undefined SCENARIO_ID');
end

% --- 3.2 CACCm gains and HDV transfer function parameters ---
if SCENARIO_ID == 1
    k1_opt = [0.08 0.08 0.19];      % proportional gain for each CAV type
    k2_opt = [1.33 1.33 0.37];      % derivative gain
    HDV_v  = [0.64 0.23 1.09;       % HDV dynamics parameters (used in feedforward)
              0.64 0.23 1.09; 
              1.02 0.7  1.14];
elseif SCENARIO_ID == 2
    k1_opt = [0.05 0.05 0.06];
    k2_opt = [0.78 0.78 1.09];
    HDV_v  = [0.62 0.33 1.14; 
              0.62 0.33 1.14; 
              0.56 0.92 1.1];
end

%% ================== 4. Precomputation and Initialization ==================
Se_HDV_vec = S_st + Ve * th;      % equilibrium spacing for HDVs [m]
Se_CAV = D + Ve * h;              % equilibrium spacing for CAVs [m]

% --- State matrices ---
Position = zeros(Num_platoon+1, total_steps);
Velocity = ones(Num_platoon+1, total_steps) * Ve;
Acceleration = zeros(Num_platoon+1, total_steps);
Jerk = zeros(Num_platoon+1, total_steps);
U = zeros(Num_platoon+1, total_steps);  % desired acceleration input

% --- Initial positions (starting from standstill) ---
for i = 1:Num_platoon
    dis_gap = 0;
    for j = 1:i
        if Vehicle_type(j) == 1
            dis_gap = dis_gap + Se_CAV;
        else
            dis_gap = dis_gap + Se_HDV_vec(j);
        end
    end
    Position(i,1) = -dis_gap;
end
Position(Num_platoon+1,1) = 0;   % leader at origin

% --- CACCm feedforward: HDV string convolution kernel ---
syms s k
F_sys = sym(zeros(1,3));
for i = 1:3
    H_v_s = (HDV_v(i,2)*s + HDV_v(i,1)/HDV_v(i,3)) / ...
            (s^2 + (HDV_v(i,1)+HDV_v(i,2))*s + HDV_v(i,1)/HDV_v(i,3));
    curr_n_hdv = Num_HDV(index == i);
    if ~isempty(curr_n_hdv)
        F_sys(i) = (tao*s+1)/(h*s+1) * H_v_s^curr_n_hdv(1);
    else
        F_sys(i) = 0;
    end
end

J_len = fix(10/deta_t) + 1;   % 10 seconds of kernel length
J = 1:J_len;
G = zeros(3, length(J));
for i = 1:3
    if F_sys(i) ~= 0
        f_time = ilaplace(F_sys(i), s, k);
        for j_idx = J
            G(i, j_idx) = double(subs(f_time, k, (j_idx-1)*deta_t));
        end
    end
end

% --- CAV-CAV feedforward filter coefficients (for leader acceleration) ---
filt_a0 = 2*h + deta_t;   % denominator coefficients
filt_a1 = deta_t - 2*h;
filt_b0 = 2*tao + deta_t; % numerator coefficients
filt_b1 = deta_t - 2*tao;
U_filter_prev = zeros(1, Num_platoon);   % previous filtered acceleration

%% ================== 5. Simulation Loop ==================
step_init_end = round(t_init / deta_t);                    % end of initialization phase
step_perturb_start = round((t_init + t_stable_run) / deta_t) + 1;  % start of perturbation

% --- Phase 1: initialization (all vehicles at equilibrium speed) ---
for t = 2 : step_init_end
    for i = 1 : Num_platoon+1
        Position(i,t) = Position(i,t-1) + Ve*deta_t;
    end
end

% --- Precompute leader trajectory for all time steps ---
for t = step_init_end + 1 : total_steps
    acc_val = 0;
    % Apply perturbation only during the defined interval
    if t >= step_perturb_start && t < step_perturb_start + length(HDV_deacceleration)
        acc_val = HDV_deacceleration(t - step_perturb_start + 1);
    end
    
    Acceleration(Num_platoon+1, t) = acc_val;
    Velocity(Num_platoon+1,t) = Velocity(Num_platoon+1,t-1) + Acceleration(Num_platoon+1,t-1)*deta_t;
    Position(Num_platoon+1,t) = Position(Num_platoon+1,t-1) + Velocity(Num_platoon+1,t-1)*deta_t;
end

% --- Car‑following control (from end of initialization onward) ---
for t = step_init_end + 1 : total_steps
    
    % --- Loop over all following vehicles ---
    for i = 1:Num_platoon
        
        % ----- 1. HDV: OVM model -----
        if Vehicle_type(i) == 0 
            % Account for driver reaction time phi(i)
            idx_prev = t - phi(i); 
            if idx_prev < 1, idx_prev = 1; end
            
            % Get preceding vehicle's state (leader if i==1)
            if i == 1
                pre_v = Velocity(Num_platoon+1, idx_prev);
                pre_p = Position(Num_platoon+1, idx_prev);
            else
                pre_v = Velocity(i-1, idx_prev);
                pre_p = Position(i-1, idx_prev);
            end
            
            delta_p = pre_p - Position(i, idx_prev);
            if delta_p > S_go
                V_opt = V_max;
            elseif delta_p < S_st
                V_opt = 0;
            else
                V_opt = (1/th(i)) * (delta_p - S_st);
            end
            
            Acceleration(i, t) = alpha(i)*(V_opt - Velocity(i, idx_prev)) + ...
                                 beta(i)*(pre_v - Velocity(i, idx_prev));
        
        % ----- 2. CAV: CACCm control -----
        else 
            % --- Gather information ---
            if i == 1
                target_front = Num_platoon+1;
            else
                target_front = i - 1;
            end
            p_front = Position(target_front, t-1);
            v_front = Velocity(target_front, t-1);
            
            % For feedforward: find the farthest vehicle whose state is used
            if i == 1
                target_far = Num_platoon+1;
            else
                target_far = i - Num_HDV(i) - 1;
                if target_far <= 0
                    target_far = Num_platoon+1;
                end
            end
            
            % --- Feedback: spacing error and its derivative ---
            e = p_front - Position(i,t-1) - h*Velocity(i,t-1) - D;
            ee = v_front - Velocity(i,t-1) - h*Acceleration(i,t-1);
            
            % Gains depending on position in platoon (index)
            if index(i) > 0
                k1_curr = k1_opt(index(i));
                k2_curr = k2_opt(index(i));
            else
                % Fallback (should not happen for CAVs)
                k1_curr = 0.5;
                k2_curr = 1.0;
            end
            
            % --- Feedforward term ---
            u_ff = 0;
            if Num_HDV(i) == 0 || i == 1
                % No HDV ahead: use leader acceleration via low‑pass filter
                acc_in_curr = Acceleration(target_far, t-theta);
                acc_in_prev = Acceleration(target_far, t-theta-1);
                u_ff = (filt_b0*acc_in_curr + filt_b1*acc_in_prev - filt_a1*U_filter_prev(i))/filt_a0;
                U_filter_prev(i) = u_ff;
            elseif Num_HDV(i) <= 3 && index(i) > 0
                % HDVs ahead: convolution with precomputed kernel
                m_len = length(Acceleration(target_far,:));
                temp = t - theta;
                J_indices = max(1, temp+1-length(J)) : min(temp, m_len);
                if temp > 0
                     conv_res = sum(Acceleration(target_far, J_indices) .* G(index(i), temp-J_indices+1));
                     u_ff = conv_res * deta_t;
                end
            end
            
            % Combined control signal
            u_control = k1_curr*e + k2_curr*ee + u_ff;
            U(i, t) = u_control;
            
            % --- Vehicle dynamics (first‑order lag) ---
            Jerk(i, t) = 1/tao * (U(i, t) - Acceleration(i, t-1));
            Acceleration(i, t) = Acceleration(i, t-1) + Jerk(i, t-1) * deta_t;
        end
        
        % --- Enforce acceleration limits ---
        if Acceleration(i, t) > max_acc
            Acceleration(i, t) = max_acc;
        elseif Acceleration(i, t) < min_acc
            Acceleration(i, t) = min_acc;
        end
    end 
    
    % --- Update velocity and position for all following vehicles ---
    for i = 1:Num_platoon
        Velocity(i,t) = Velocity(i,t-1) + Acceleration(i,t-1) * deta_t;
        Position(i,t) = Position(i,t-1) + Velocity(i,t-1) * deta_t;
    end
    
end 

%% ================== 6. Results and Visualization ==================
Velocity = real(Velocity); 
Acceleration = real(Acceleration);
Position = real(Position);

% --- Time window for metric computation (from perturbation start to end) ---
calc_indices = step_perturb_start : total_steps;
num_calc_steps = length(calc_indices);
duration_calc = num_calc_steps * deta_t;

Vel_calc = Velocity(:, calc_indices);
Acc_calc = Acceleration(:, calc_indices);
Pos_calc = Position(:, calc_indices);

% --- Initialize metric arrays ---
PAV_vec = zeros(Num_platoon, 1);   % Peak Acceleration Variation
AAV_vec = zeros(Num_platoon, 1);   % Average Acceleration Variation
PSV_vec = zeros(Num_platoon, 1);   % Peak Speed Variation
DRAC_vec = zeros(Num_platoon, 1);  % Deceleration Rate to Avoid Crash
Fuel_vec = zeros(Num_platoon, 1);  % Fuel Consumption Rate (FCR)

for num = 1:Num_platoon
    v_ego = Vel_calc(num, :);
    a_ego = Acc_calc(num, :);
    p_ego = Pos_calc(num, :);
    
    if num == 1
        v_front = Vel_calc(Num_platoon+1, :);
        p_front = Pos_calc(Num_platoon+1, :);
    else
        v_front = Vel_calc(num-1, :);
        p_front = Pos_calc(num-1, :);
    end
    
    % --- (1) PAV ---
    PAV_vec(num) = max(abs(a_ego));
    
    % --- (2) AAV ---
    AAV_vec(num) = sum(abs(a_ego)) / num_calc_steps;
    
    % --- (3) PSV (deviation from equilibrium speed) ---
    PSV_vec(num) = max(abs(v_ego - Ve));
    
    % --- (4) DRAC ---
    rel_v = v_ego - v_front;   % relative speed (ego – front)
    gap = p_front - p_ego;
    gap(gap < 0.1) = 0.1;      % avoid division by zero
    drac_vals = (rel_v.^2) ./ gap;
    drac_vals(rel_v <= 0) = 0; % only when closing in
    DRAC_vec(num) = max(drac_vals);
    
    % --- (5) Fuel Consumption Rate (based on Virginia Tech model) ---
    ff = zeros(1, num_calc_steps);
    for k = 1:num_calc_steps
        v_k = v_ego(k);
        a_k = a_ego(k);
        RR = 0.333 + 0.00108*v_k^2 + 1.2*a_k;
        if RR > 0
            if a_k > 0
                ff(k) = 0.444 + 0.09*RR*v_k + 0.054*a_k^2*v_k;
            else
                ff(k) = 0.444 + 0.09*RR*v_k;
            end
        else
            ff(k) = 0.444;
        end
    end
    Fuel_vec(num) = sum(ff) * deta_t / duration_calc;   % [mL/s]
end

% --- Mean values over all vehicles and over CAVs only ---
Mean_All = [mean(PAV_vec), mean(AAV_vec), mean(PSV_vec), mean(DRAC_vec), mean(Fuel_vec)];

cav_indices = find(Vehicle_type == 1);
if ~isempty(cav_indices)
    Mean_CAV = [mean(PAV_vec(cav_indices)), mean(AAV_vec(cav_indices)), ...
                mean(PSV_vec(cav_indices)), mean(DRAC_vec(cav_indices)), ...
                mean(Fuel_vec(cav_indices))];
else
    Mean_CAV = zeros(1,5);
end

% --- Display results ---
disp('------------------------------------------------------------');
disp(['Strategy: CACCm | Scenario: Type ' num2str(SCENARIO_ID)]);
disp('Time window: from perturbation start to simulation end');
disp('Metrics: [MAV, AAV, MSV, DRAC, FCR]');
disp(' ');
disp('1. Average over all vehicles (10):');
disp(Mean_All);
disp(' ');
disp(['2. Average over CAVs only (' num2str(length(cav_indices)) ' vehicles):']);
disp(Mean_CAV);
disp('------------------------------------------------------------');

% --- Plotting ---
figure('Name', ['CACCm - Type ' num2str(SCENARIO_ID)], 'Color', 'w', 'Position', [100, 100, 1000, 700]);
t_axis = (1:total_steps) * deta_t;
Color_HDV = [1 0.8 0.2];
Color_CAV = [0 0.2 0.8];

% --- Subplot 1: full simulation ---
subplot(2, 1, 1);
hold on; box on; grid on;
plot(t_axis, Velocity(Num_platoon+1,:), 'k', 'LineWidth', 2, 'DisplayName', 'Leader');
for i = 1:Num_platoon
    if Vehicle_type(i) == 1
        c = Color_CAV; stl = '-';
    else
        c = Color_HDV; stl = '-.';
    end
    if i == find(Vehicle_type==1, 1) || i == find(Vehicle_type==0, 1)
        plot(t_axis, Velocity(i,:), 'Color', c, 'LineStyle', stl, 'LineWidth', 1.5, ...
             'DisplayName', ['Type ' num2str(Vehicle_type(i))]);
    else
        plot(t_axis, Velocity(i,:), 'Color', c, 'LineStyle', stl, 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    end
end
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title(['CACCm - Full Simulation (Type ' num2str(SCENARIO_ID) ')']);
legend('Leader', 'CAV', 'HDV', 'Location', 'best');
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

% --- Subplot 2: zoom on perturbation phase ---
subplot(2, 1, 2);
hold on; box on; grid on;
start_plot_idx = step_perturb_start;
t_sub = t_axis(start_plot_idx:end);
v_sub = Velocity(:, start_plot_idx:end);

plot(t_sub, v_sub(Num_platoon+1,:), 'k', 'LineWidth', 2, 'DisplayName', 'Leader');
for i = 1:Num_platoon
    if Vehicle_type(i) == 1
        c = Color_CAV; stl = '-';
    else
        c = Color_HDV; stl = '-.';
    end
    if i == find(Vehicle_type==1, 1) || i == find(Vehicle_type==0, 1)
        plot(t_sub, v_sub(i,:), 'Color', c, 'LineStyle', stl, 'LineWidth', 1.5, ...
             'DisplayName', ['Type ' num2str(Vehicle_type(i))]);
    else
        plot(t_sub, v_sub(i,:), 'Color', c, 'LineStyle', stl, 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    end
end
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title(['Perturbation Phase (from ' num2str(t_stable_run + t_init) ' s onward)']);
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
xlim([t_sub(1), t_sub(end)]);