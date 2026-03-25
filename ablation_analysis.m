clear;
clc;

%% ================= 1. Configuration =================
leader_idx = 4; % Leader index
target_idx = 3; % Target index
m_hdv = leader_idx - target_idx; 
if m_hdv < 1, error('Leader index must be greater than Target index'); end

Num_Subject_CAV = 1;
deta_t = 0.1;

Method_List = {'ACC (FB Only)', 'Fixed CACCm', 'Optimized CACCm'};
Num_Methods = length(Method_List);

%% ================= 2. Data Loading =================
try
    load acc_ngsim.mat
    load vel_ngsim.mat
    load pos_ngsim.mat

    start_step = 1;
    stop_step = size(acceleration, 2);

    acc_target_real = acceleration(target_idx, start_step:stop_step);
    acc_leader_real = acceleration(leader_idx, start_step:stop_step);
    vel_leader_real = velocity(leader_idx, start_step:stop_step);
    pos_leader_real = position(leader_idx, start_step:stop_step);
    
    total_steps = length(acc_leader_real);
    T_sim = total_steps * deta_t;
catch
    error('Data files not found.');
end

%% ================= 3. Initialization & Storage =================
h = 1.2; D = 5; tao = 0.5; theta = round(0.1/deta_t); 
V_max = 30; max_acc = 5; min_acc = -5;

Result_Matrix = zeros(Num_Methods, 5); 

Velocity_Store = zeros(Num_Methods, total_steps); 

%% ================= 4. Ablation Loop =================
fprintf('Start Ablation Study for m=%d, n=%d...\n', m_hdv, Num_Subject_CAV);

for m_idx = 1:Num_Methods
    method_name = Method_List{m_idx};
    fprintf('  Running Method %d: %s ... ', m_idx, method_name);
    
    % --- 4.1 parameters (The Switchboard) ---
    k1_curr = 0; k2_curr = 0; hv_curr = [0 0 0];
    
    switch m_idx
        case 1 % === ACC (Feedback Only) ===
            k1_curr = 0.5; 
            k2_curr = 1.0;
            hv_curr = [0, 0, 0];
            
        case 2 % === Fixed CACCm (No Optimization) ===
            k1_curr = 0.5; 
            k2_curr = 1.0;
            hv_curr = [0.2, 0.4, 1.5];
            
        case 3 % === Optimized CACCm (Full Method) ===
            if m_hdv == 1
                if Num_Subject_CAV == 1
                    hv_curr = [0.64, 0.71, 1.09]; k1_curr = 0.05; k2_curr = 1.35;
                elseif Num_Subject_CAV == 2
                    hv_curr = [0.62, 0.33, 1.14]; k1_curr = 0.05; k2_curr = 0.78;
                elseif Num_Subject_CAV >= 3
                    hv_curr = [0.56, 0.92, 1.10]; k1_curr = 0.06; k2_curr = 1.09; 
                end
            elseif m_hdv == 2
                if Num_Subject_CAV == 1
                    hv_curr = [0.64, 0.23, 1.09]; k1_curr = 0.08; k2_curr = 1.33;
                elseif Num_Subject_CAV >= 2
                    hv_curr = [0.63, 0.47, 1.15]; k1_curr = 0.02; k2_curr = 1.81; 
                end
            elseif m_hdv == 3
                if Num_Subject_CAV == 1
                    hv_curr = [1.02, 0.70, 1.14]; k1_curr = 0.19; k2_curr = 0.37; 
                else
                    % Fallback
                    hv_curr = [1.02, 0.70, 1.14]; k1_curr = 0.5; k2_curr = 1;
                end
            else
                % Fallback defaults
                hv_curr = [0.2, 0.4, 1.5]; k1_curr = 0.5; k2_curr = 1.0;
            end
    end
    
    % --- 4.2  G(t) ---
    G_kernel = [];
    if sum(hv_curr) > 0
        syms s k_sym
        alpha_v = hv_curr(1); beta_v = hv_curr(2); th_v = hv_curr(3);
        H_v_s = (beta_v*s + alpha_v/th_v) / (s^2 + (alpha_v + beta_v)*s + alpha_v/th_v);
        F_s = ((tao*s + 1) / (h*s + 1)) * (H_v_s)^m_hdv;
        f_t = ilaplace(F_s, s, k_sym);
        
        J_len = fix(10/deta_t) + 1; 
        G_kernel = zeros(1, J_len);
        for j = 1:J_len
            t_val = (j-1)*deta_t;
            G_kernel(j) = double(subs(f_t, k_sym, t_val));
        end
    end
    
    % --- 4.3 simulation ---
    Pos_CAV = zeros(Num_Subject_CAV, total_steps);
    Vel_CAV = zeros(Num_Subject_CAV, total_steps);
    Acc_CAV = zeros(Num_Subject_CAV, total_steps);
    U_CAV = zeros(Num_Subject_CAV, total_steps);

    for i = 1:Num_Subject_CAV
        if i == 1
            Pos_CAV(i,1) = pos_leader_real(1) - (D + vel_leader_real(1)*h);
            Vel_CAV(i,1) = vel_leader_real(1);
        else
            Pos_CAV(i,1) = Pos_CAV(i-1,1) - (D + Vel_CAV(i-1,1)*h);
            Vel_CAV(i,1) = Vel_CAV(i-1,1);
        end
    end
    
    for t = 2 : total_steps
        for i = 1 : Num_Subject_CAV
            if i == 1
                p_front = pos_leader_real(t-1); v_front = vel_leader_real(t-1); a_front = acc_leader_real(t-1);
            else
                p_front = Pos_CAV(i-1, t-1); v_front = Vel_CAV(i-1, t-1); a_front = Acc_CAV(i-1, t-1);
            end
            p_ego = Pos_CAV(i, t-1); v_ego = Vel_CAV(i, t-1); a_ego = Acc_CAV(i, t-1);
            
            
            if i == 1
                e = p_front - p_ego - h*v_ego - D; 
                ee = v_front - v_ego - h*a_ego;
                
                u_ff = 0;
                if ~isempty(G_kernel)
                    t_eff = t - theta;
                    idx_end = t_eff;
                    idx_start = t_eff - length(G_kernel) + 1;
                    range_indices = idx_start : idx_end;
                    valid_mask = range_indices > 0;
                    valid_indices = range_indices(valid_mask);
                    if ~isempty(valid_indices)
                        acc_chunk = acc_target_real(valid_indices);
                        kernel_indices = t_eff - valid_indices + 1;
                        g_chunk = G_kernel(kernel_indices);
                        u_ff = sum(acc_chunk .* g_chunk) * deta_t;
                    end
                end
                
                u_val = k1_curr*e + k2_curr*ee + u_ff;
                
            else
                e = p_front - p_ego - h*v_ego - D; 
                ee = v_front - v_ego - h*a_ego;
                
                if m_idx == 1 
                    u_val = 0.5*e + 1.0*ee; 
                else
                    u_val = 0.5*e + 1.0*ee + a_front; 
                end
            end
            
            Jerk = (1/tao) * (u_val - a_ego);
            Acc_CAV(i,t) = a_ego + Jerk * deta_t;
            if Acc_CAV(i,t) > max_acc, Acc_CAV(i,t) = max_acc; end
            if Acc_CAV(i,t) < min_acc, Acc_CAV(i,t) = min_acc; end
            Vel_CAV(i,t) = max(0, v_ego + Acc_CAV(i,t-1) * deta_t);
            Pos_CAV(i,t) = p_ego + Vel_CAV(i,t-1) * deta_t;
        end
    end
    
    Vel_CAV = real(Vel_CAV); Acc_CAV = real(Acc_CAV); Pos_CAV = real(Pos_CAV);
    
    % --- 4.4 Evaluation ---
    idx_eval = Num_Subject_CAV;
    
    Velocity_Store(m_idx, :) = Vel_CAV(idx_eval, :);
    
    % 1. MAV
    Result_Matrix(m_idx, 1) = max(abs(Acc_CAV(idx_eval,:))); 
    % 2. AAV
    Result_Matrix(m_idx, 2) = sum(abs(Acc_CAV(idx_eval,:))) / total_steps; 
    % 3. MSV
    Result_Matrix(m_idx, 3) = max(abs(Vel_CAV(idx_eval,:) - mean(Vel_CAV(idx_eval,:)))); 
    % 4. DRAC
    if idx_eval == 1
        v_pre = vel_leader_real; p_pre = pos_leader_real;
    else
        v_pre = Vel_CAV(idx_eval-1,:); p_pre = Pos_CAV(idx_eval-1,:);
    end
    rel_v = Vel_CAV(idx_eval,:) - v_pre;
    gap = p_pre - Pos_CAV(idx_eval,:);
    gap(gap<0.1) = 0.1;
    drac_vals = (rel_v.^2) ./ gap;
    drac_vals(rel_v <= 0) = 0;
    Result_Matrix(m_idx, 4) = max(drac_vals); 
    
    % 5. FCR
    ff = zeros(1, total_steps);
    for k=1:total_steps
        v_k = Vel_CAV(idx_eval,k); a_k = Acc_CAV(idx_eval,k);
        RR = 0.333 + 0.00108*v_k^2 + 1.2*a_k;
        if RR>0
            if a_k>0, ff(k) = 0.444 + 0.09*RR*v_k + 0.054*a_k^2*v_k;
            else,     ff(k) = 0.444 + 0.09*RR*v_k; end
        else, ff(k) = 0.444; end
    end
    Result_Matrix(m_idx, 5) = sum(ff)*deta_t / T_sim;
    
    fprintf('Done.\n');
end

%% ================= 5. Results & Visualization =================
disp('======================================================');
disp(['Ablation Results (Last CAV, n=' num2str(Num_Subject_CAV) ')']);
disp('cols: [MAV, AAV, MSV, DRAC, FCR]');
disp('------------------------------------------------------');
disp('Method 1: ACC (FB Only, Fixed Gains)');
disp(Result_Matrix(1,:));
disp('Method 2: ACC + Feedforward (Fixed Params)');
disp(Result_Matrix(2,:));
disp('Method 3: CACCm (Optimized Params - Proposed)');
disp(Result_Matrix(3,:));
disp('======================================================');

% --- Plotting Trajectories ---
t_axis = (0:total_steps-1) * deta_t;
figure('Name', 'Ablation Study Trajectories', 'Color', 'w');
hold on; box on; grid on;

% 1. Plot Leader
p_lead = plot(t_axis, vel_leader_real, 'k--', 'LineWidth', 1.5);

% 2. Plot Methods
colors = lines(3); 
% Method 1: ACC (Blueish)
p1 = plot(t_axis, Velocity_Store(1,:), 'Color', [0 0.447 0.741], 'LineWidth', 1.5, 'LineStyle', ':');
% Method 2: Fixed CACC (Reddish/Orange)
p2 = plot(t_axis, Velocity_Store(2,:), 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5, 'LineStyle', '-.');
% Method 3: Optimized CACCm (Greenish - Proposed)
p3 = plot(t_axis, Velocity_Store(3,:), 'Color', [0.466 0.674 0.188], 'LineWidth', 2.0, 'LineStyle', '-');

xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontName', 'Times New Roman');
title(['Ablation Study: Velocity Profile of Last CAV (n=' num2str(Num_Subject_CAV) ')'], 'FontSize', 12);

legend([p_lead, p1, p2, p3], ...
    {'Leader', 'ACC (FB Only)', 'Fixed CACCm', 'Optimized CACCm (Proposed)'}, ...
    'Location', 'Best', 'FontSize', 10);