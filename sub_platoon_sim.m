clear;
clc;

%% ================= 1. Configuration =================
% Core settings
leader_idx = 4;      % Leader (NGSIM, immediate predecessor)
target_idx = 3;      % Target (NGSIM, far predecessor)

% Calculate m (number of HDVs crossed)
m_hdv = leader_idx - target_idx;
if m_hdv < 1
    error('Leader index must be greater than Target index');
end

% Simulation parameters
Num_Subject_CAV = 1;   % n = 1 or 2
deta_t = 0.1;          % Time step [s]

%% ================= 2. Data Loading =================
% Load NGSIM trajectories
try
    load acc_ngsim.mat
    load vel_ngsim.mat
    load pos_ngsim.mat

    start_step = 1;
    stop_step = size(acceleration, 2);

    % Extract real trajectories
    acc_target_real = acceleration(target_idx, start_step:stop_step);
    vel_target_real = velocity(target_idx, start_step:stop_step);

    acc_leader_real = acceleration(leader_idx, start_step:stop_step);
    vel_leader_real = velocity(leader_idx, start_step:stop_step);
    pos_leader_real = position(leader_idx, start_step:stop_step);

    vel_real_all = velocity(:, start_step:stop_step);

    total_steps = length(acc_leader_real);
    T_sim = total_steps * deta_t;

catch
    error('Data files not found.');
end

%% ================= 3. Parameter Initialization =================
% General parameters
h = 1.2;          % Time gap [s]
D = 5;            % Standstill distance [m]
tao = 0.5;        % Actuator lag [s]
theta = round(0.1 / deta_t);  % Communication delay [steps]
V_max = 30;       % Maximum velocity [m/s]
max_acc = 5;      % Maximum acceleration [m/s²]
min_acc = -5;     % Minimum acceleration [m/s²]

% CACCm gains (depending on m and n)
k1_opt = 0.5; k2_opt = 1.0; hv_params = [0 0 0];

% Select parameters based on m_hdv and Num_Subject_CAV
if m_hdv == 1
    if Num_Subject_CAV == 1
        hv_params = [0.64, 0.71, 1.09];
        k1_opt = 0.05; k2_opt = 1.35;
    elseif Num_Subject_CAV == 2
        hv_params = [0.62, 0.33, 1.14];
        k1_opt = 0.05; k2_opt = 0.78;
    elseif Num_Subject_CAV >= 3
        hv_params = [0.56, 0.92, 1.10];
        k1_opt = 0.06; k2_opt = 1.09;
    end
elseif m_hdv == 2
    if Num_Subject_CAV == 1
        hv_params = [0.64, 0.23, 1.09];
        k1_opt = 0.08; k2_opt = 1.33;
    elseif Num_Subject_CAV >= 2
        hv_params = [0.63, 0.47, 1.15];
        k1_opt = 0.02; k2_opt = 1.81;
    end
elseif m_hdv == 3
    if Num_Subject_CAV == 1
        hv_params = [1.02, 0.70, 1.14];
        k1_opt = 0.19; k2_opt = 0.37;
    end
end

%% ================= 4. Feedforward Kernel Calculation =================
fprintf('Calculating feedforward kernel for CACCm (m=%d, n=%d)...\n', m_hdv, Num_Subject_CAV);
syms s k_sym
alpha_v = hv_params(1);
beta_v  = hv_params(2);
th_v    = hv_params(3);

% Virtual HDV transfer function H(s)
H_v_s = (beta_v * s + alpha_v / th_v) / (s^2 + (alpha_v + beta_v) * s + alpha_v / th_v);

% Feedforward filter: G_cav(s) = (tao*s+1)/(h*s+1)
G_cav_s = (tao * s + 1) / (h * s + 1);

% Overall feedforward: F(s) = G_cav(s) * (H_v(s))^m
F_s = G_cav_s * (H_v_s)^m_hdv;

% Inverse Laplace transform to obtain impulse response f(t)
f_t = ilaplace(F_s, s, k_sym);

% Discretize to obtain kernel G(k)
J_len = fix(10 / deta_t) + 1;
G_kernel = zeros(1, J_len);
for j = 1:J_len
    t_val = (j-1) * deta_t;
    G_kernel(j) = double(subs(f_t, k_sym, t_val));
end
fprintf('Kernel calculation complete.\n');

%% ================= 5. State Initialization =================
Pos_CAV = zeros(Num_Subject_CAV, total_steps);
Vel_CAV = zeros(Num_Subject_CAV, total_steps);
Acc_CAV = zeros(Num_Subject_CAV, total_steps);
Jerk_CAV = zeros(Num_Subject_CAV, total_steps);
U_CAV = zeros(Num_Subject_CAV, total_steps);   % Desired acceleration

% Initial positions and velocities (assuming constant spacing)
for i = 1:Num_Subject_CAV
    if i == 1
        Pos_CAV(i,1) = pos_leader_real(1) - (D + vel_leader_real(1) * h);
        Vel_CAV(i,1) = vel_leader_real(1);
    else
        Pos_CAV(i,1) = Pos_CAV(i-1,1) - (D + Vel_CAV(i-1,1) * h);
        Vel_CAV(i,1) = Vel_CAV(i-1,1);
    end
end

%% ================= 6. Simulation Loop =================
fprintf('Starting simulation: CACCm | m=%d | n=%d\n', m_hdv, Num_Subject_CAV);

for t = 2:total_steps
    for i = 1:Num_Subject_CAV

        % ---- A. Perception ----
        if i == 1
            p_front = pos_leader_real(t-1);
            v_front = vel_leader_real(t-1);
            a_front = acc_leader_real(t-1);
        else
            p_front = Pos_CAV(i-1, t-1);
            v_front = Vel_CAV(i-1, t-1);
            a_front = Acc_CAV(i-1, t-1);
        end

        % Far predecessor information (only for the first CAV)
        if t > theta
            v_far = vel_target_real(t - theta);
            a_far = acc_target_real(t - theta);
        else
            v_far = vel_target_real(1);
            a_far = acc_target_real(1);
        end

        p_ego = Pos_CAV(i, t-1);
        v_ego = Vel_CAV(i, t-1);
        a_ego = Acc_CAV(i, t-1);

        % ---- B. Control Law ----
        if i == 1
            % First CAV uses the full CACCm strategy
            % Spacing error and relative velocity error
            e  = p_front - p_ego - h * v_ego - D;
            ee = v_front - v_ego - h * a_ego;

            % Feedforward term: convolution of far predecessor acceleration with kernel
            u_ff = 0;
            if ~isempty(G_kernel)
                t_eff = t - theta;                     % Effective time index
                idx_end = t_eff;
                idx_start = t_eff - length(G_kernel) + 1;
                range_indices = idx_start:idx_end;
                valid_mask = range_indices > 0;
                valid_indices = range_indices(valid_mask);

                if ~isempty(valid_indices)
                    acc_chunk = acc_target_real(valid_indices);
                    kernel_indices = t_eff - valid_indices + 1;
                    g_chunk = G_kernel(kernel_indices);
                    u_ff = sum(acc_chunk .* g_chunk) * deta_t;
                end
            end

            u_val = k1_opt * e + k2_opt * ee + u_ff;

        else
            % Following CAVs use standard CACC with predecessor feedforward
            e  = p_front - p_ego - h * v_ego - D;
            ee = v_front - v_ego - h * a_ego;
            u_val = 0.5 * e + 1.0 * ee + a_front;
        end

        % ---- C. Vehicle Dynamics (First‑order lag) ----
        U_CAV(i, t) = u_val;
        Jerk_CAV(i, t) = (1 / tao) * (U_CAV(i, t) - Acc_CAV(i, t-1));
        Acc_CAV(i, t) = Acc_CAV(i, t-1) + Jerk_CAV(i, t) * deta_t;

        % Limit acceleration
        if Acc_CAV(i, t) > max_acc
            Acc_CAV(i, t) = max_acc;
        end
        if Acc_CAV(i, t) < min_acc
            Acc_CAV(i, t) = min_acc;
        end

        % Update velocity and position (Euler integration)
        Vel_CAV(i, t) = Vel_CAV(i, t-1) + Acc_CAV(i, t-1) * deta_t;
        if Vel_CAV(i, t) < 0
            Vel_CAV(i, t) = 0;
        end
        Pos_CAV(i, t) = Pos_CAV(i, t-1) + Vel_CAV(i, t-1) * deta_t;
    end
end

% Ensure real parts only (symbolic computation may produce small imaginary parts)
Vel_CAV = real(Vel_CAV);
Acc_CAV = real(Acc_CAV);
Pos_CAV = real(Pos_CAV);

%% ================= 7. Performance Metrics =================
idx_eval = Num_Subject_CAV;   % Evaluate the last CAV in the platoon
metrics = zeros(1, 5);

% Peak acceleration variation (PAV)
metrics(1) = max(abs(Acc_CAV(idx_eval, :)));

% Average absolute acceleration variation (AAV)
metrics(2) = sum(abs(Acc_CAV(idx_eval, :))) / total_steps;

% Peak speed variation (PSV)
metrics(3) = max(abs(Vel_CAV(idx_eval, :) - mean(Vel_CAV(idx_eval, :))));

% Deceleration rate to avoid crash (DRAC)
if idx_eval == 1
    v_pre = vel_leader_real;
    p_pre = pos_leader_real;
else
    v_pre = Vel_CAV(idx_eval-1, :);
    p_pre = Pos_CAV(idx_eval-1, :);
end
rel_v = Vel_CAV(idx_eval, :) - v_pre;
gap   = p_pre - Pos_CAV(idx_eval, :);
gap(gap < 0.1) = 0.1;          % Avoid division by zero
metrics(4) = max(rel_v.^2 ./ gap);

% Fuel consumption rate (FCR)
ff = zeros(1, total_steps);
for k = 1:total_steps
    v_k = Vel_CAV(idx_eval, k);
    a_k = Acc_CAV(idx_eval, k);
    RR = 0.333 + 0.00108 * v_k^2 + 1.2 * a_k;
    if RR > 0
        if a_k > 0
            ff(k) = 0.444 + 0.09 * RR * v_k + 0.054 * a_k^2 * v_k;
        else
            ff(k) = 0.444 + 0.09 * RR * v_k;
        end
    else
        ff(k) = 0.444;
    end
end
metrics(5) = sum(ff) * deta_t / T_sim;   % Average fuel consumption rate

% Display results
disp('------------------------------------------------');
fprintf('Strategy: CACCm (m=%d, n=%d)\n', m_hdv, Num_Subject_CAV);
disp('Evaluated vehicle: Last CAV');
disp('Metrics:    MAV      AAV      MSV      DRAC     FCR');
fprintf('Values: %8.4f %8.4f %8.4f %8.4f %8.4f\n', metrics);
disp('------------------------------------------------');

%% ================= 8. Visualization =================
t_axis = (0:total_steps-1) * deta_t;
figure('Name', 'CACCm Simulation Results', 'Color', 'w');

plot(t_axis, Vel_CAV(end,:), 'b', 'LineWidth', 1.5);
hold on;
plot(t_axis, vel_leader_real, 'k--', 'LineWidth', 1);
ylabel('Velocity [m/s]');
xlabel('Time [s]');
title('Last CAV vs. Leader');
legend('CAV', 'Leader', 'Location', 'best');
grid on;

fprintf('Simulation finished.\n');