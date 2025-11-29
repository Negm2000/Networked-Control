clc
clear
close all

% Project 7: Coupled Penduli
% Reference: 07_Penduli.pdf
% Two inverted penduli connected by a spring (stiffness k)
k = 0.2;   % N/m (Spring stiffness)
% k = 2;   % N/m
% k = 200; % N/m
l = 1;   % m (Rod length)
m = 1;   % kg (Mass)
g = 9.8; % m/s2 (Gravity)
a = l;   % Spring attachment height
h = 0.5; % s (Sampling time)
N = 2;   % Number of subsystems

% Continuous-Time State Space Matrices (A, B, C)
% State x = [theta1, d_theta1, theta2, d_theta2]'
A = [
    0               1       0         0
    g/l-k*a^2/m/l^2 0 k*a^2/m/l^2     0
    0               0      0          1
    k*a^2/m/l^2     0 g/l-k*a^2/m/l^2 0
    ];
B=[[0 1/m/l^2 0 0]',[0 0 0 1/m/l^2]'];
C=eye(4); % Full state feedback assumed

%  1. Modelling & Discretization
% ---------------------------------
% Decompose the state and input vectors into subvectors
% Subsystem 1: States 1,2; Input 1
% Subsystem 2: States 3,4; Input 2
B1 = B(:,1);
B2 = B(:,2);
C1 = C(1:2,:);
C2 = C(3:end,:);

% Discretization 
% F = e^(Ah)
F = expm(A*h); 
% G = integral(e^(At)dt) * B
G = A \ (expm(A*h) - expm(A*0)) * B; 
H = C;

% Decomposed discrete matrices
G1 = G(:,1);
G2 = G(:,2);
H1 = H(1:2,:);
H2 = H(3:end,:);

% 2. Analysis
% ---------------------------------
% Check open-loop stability
% Unstable if any eigenvalue magnitude > 1
isStable(F);

% 3. Design - LMI Setup
n_states = 4;
m_inputs = 2;

% Optimization Variables
% P: Lyapunov matrix (P > 0)
% L: Controller proxy variable (L = K*P) to linearize the inequality
P=sdpvar(n_states);
L=sdpvar(m_inputs,n_states);

%  A. Static State-Feedback Stability
% ---------------------------------
% Condition: (F+GK)P(F+GK)' - P < 0
% Implemented using Schur Complement:
LMIconstr_stability=[[P-F*P*F'-F*L'*G'-G*L*F' , G*L;
                                L'*G'         , P  ] >= 1e-2*eye(n_states*2)];

Kx = solveLMI(LMIconstr_stability,P,L);
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable (Stability LMI).")
isStable(F_cl);

% VISUALIZATION A
generate_plots_smooth(Kx, A, B, h, n_states, 'Stability LMI');

%  B. Pole Placement (Spectral Radius / Speed)
% ---------------------------------
% Goal: Force eigenvalues |lambda| < rho
% Require (1/rho)*(F+GK) to be Schur stable.
% This scales the P term by rho^2 in the LMI.
rho = 0.12; 
LMIconstr_speed=[[rho^2*P-F*P*F'-F*L'*G'-G*L*F' , G*L;
                                L'*G'           , P  ] >= 1e-2*eye(n_states*2)];

Kx = solveLMI(LMIconstr_speed,P,L);
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable (Speed LMI).")
isStable(F_cl);

% VISUALIZATION B
generate_plots_smooth(Kx, A, B, h, n_states, 'Speed LMI');

%  C. H2 Optimal Design 
% ---------------------------------
% Goal: Minimize Trace(S), where S is an upper bound on output energy z_k.
P=sdpvar(n_states); 
L=sdpvar(m_inputs,n_states);

% Bryson's Rule for Weights
Q = diag([0.01^-2, 0.1^-2, 0.01^-2, 0.1^-2]); % State penalty
R = 10*eye(2); % Input penalty

% Generalized Plant Matrices
Gw = eye(4); % Noise input matrix (assume noise on all states)
H = [sqrtm(Q); zeros(2,4)]; % Performance State weight
Du = [zeros(4,2); sqrtm(R)]; % Performance Input weight

S = sdpvar(6, 6, 'symmetric'); 

% LMI 1: Reachability Gramian Constraint (Stability with noise)
% Ensures state energy P is bounded in presence of noise Gw
LMI_1 = [ [P - F*P*F' - F*L'*G' - G*L*F' - Gw*Gw',  G*L];
          [L'*G',                                   P] ] >= 1e-6 * eye(n_states*2);

% LMI 2: Observability/Cost Constraint
% Connects state energy P to cost S using Schur Complement
LMI_2 = [ [S,                           H*P + Du*L];
          [(H*P + Du*L)',     P] ] >= 1e-6 * eye(2*n_states+m_inputs);

LMIconstr_H2 = [LMI_1, LMI_2];
options = sdpsettings('verbose', 0); 

% Minimize Trace(S) ==> Minimize H2 Norm
J = optimize(LMIconstr_H2,trace(S),options);

if J.problem 
    disp("Unfeasible")
end

L_val = double(L);
P_val = double(P);
Kx = L_val/P_val;

F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable (H2).")
isStable(F_cl);

% VISUALIZATION C
generate_plots_smooth(Kx, A, B, h, n_states, 'H2 Design');

%  D. H-Infinity Robust Design
% ---------------------------------
% Goal: Minimize worst-case gain gamma from disturbance w to output z.
% Uses the Bounded Real Lemma in LMI form.

P = sdpvar(n_states);
L = sdpvar(m_inputs,n_states);
gamma = sdpvar(1); % Performance level (scalar)

H_perf = H;
Du_perf = Du;
Dw = zeros(size(H_perf,1), size(Gw,2)); % No direct feedthrough

% Build the 4x4 Block Matrix for Bounded Real Lemma (Page 86)
LMI_11 = P;
LMI_12 = F*P + G*L;
LMI_13 = Gw;
LMI_14 = zeros(n_states, size(H_perf,1));

LMI_22 = P;
LMI_23 = zeros(n_states, size(Gw,2));
LMI_24 = (H_perf*P + Du_perf*L)'; 

LMI_33 = gamma * eye(size(Gw,2));
LMI_34 = Dw'; 

LMI_44 = gamma * eye(size(H_perf,1));

LMI_Hinf = [ LMI_11, LMI_12, LMI_13, LMI_14;
             LMI_12',LMI_22, LMI_23, LMI_24;
             LMI_13',LMI_23',LMI_33, LMI_34;
             LMI_14',LMI_24',LMI_34',LMI_44 ];
             
constraints = [LMI_Hinf >= 1e-6 * eye(size(LMI_Hinf,1)), P >= 1e-6 * eye(n_states)];

disp('Solving H-Infinity LMI... (This may take a moment)');
sol = optimize(constraints, gamma, options); % Minimize gamma

if sol.problem == 0
    K_hinf = double(L) / double(P);
    disp('H-Infinity Controller (K_hinf) successfully designed.');
    fprintf('Optimal H-Infinity Gain (gamma): %.4f\n', double(gamma));
    
    F_cl_hinf = F + G * K_hinf;
    rho_hinf = max(abs(eig(F_cl_hinf)));
    fprintf('Closed-loop spectral radius: %.4f\n', rho_hinf);
    
    % VISUALIZATION D
    generate_plots_smooth(K_hinf, A, B, h, n_states, 'H-Inf Design');
else
    disp('H-Infinity LMI was Infeasible');
    disp(sol.info);
end

% =========================================================================
% FUNCTIONS 
% =========================================================================

function Kx = solveLMI(LMIconstr,P,L)
    % Helper function to solve standard feasibility LMI
    options = sdpsettings('verbose', 0); 
    J = optimize(LMIconstr,[],options);
    if J.problem 
        disp("Unfeasible")
        Kx = [];
        return 
    end
    L = double(L);
    P = double(P);
    Kx = L/P;
end

function rho = isStable(F)
    % Checks Spectral Radius (Theorem 3, Part I)
    rho = max(abs(eig(F)));
    if (rho > 1) 
        fprintf ('The absolute spectral radius is %.2f > 1 (Unstable).\n',rho)
    else
        fprintf("The system is stable, spectral radius = |%.2f| < 1\n", rho)
    end
end

% --- HYBRID SMOOTH PLOTTING ---
function generate_plots_smooth(K, A, B, h, n_states, plotTitle)
    % Simulates Hybrid System: 
    % - Continuous Plant (ODE)
    % - Discrete Controller (Sample & Hold)
    
    dt_physics = 0.01; % Fine integration step
    n_steps = 20; 
    t_final = h * n_steps;
    t_smooth = 0:dt_physics:t_final;
    
    x0 = [pi/6; 0; -pi/6; 0]; 
    current_x = x0;
    
    x_hist_smooth = zeros(n_states, length(t_smooth));
    u_hist = zeros(2, length(t_smooth));
    
    u_current = [0;0];
    next_sample_time = 0;
    
    for i = 1:length(t_smooth)
        t_now = t_smooth(i);
        
        % Discrete Control Update (Zero Order Hold)
        if t_now >= next_sample_time
            if ~isempty(K)
                u_current = K * current_x; 
            else
                u_current = [0;0];
            end
            next_sample_time = next_sample_time + h;
        end
        
        % Continuous Physics Update (Euler integration for simplicity)
        dx = A * current_x + B * u_current;
        current_x = current_x + dx * dt_physics;
        
        x_hist_smooth(:, i) = current_x;
        u_hist(:, i) = u_current;
    end
    
    figure('Name', plotTitle, 'Color', 'w');
    sgtitle([plotTitle ' (Smooth Simulation)']); 
    
    titles = {'Pos 1 (x1)', 'Vel 1 (x2)', 'Pos 2 (x3)', 'Vel 2 (x4)'};
    for k = 1:4
        subplot(2, 2, k);
        plot(t_smooth, x_hist_smooth(k, :), 'LineWidth', 1.5);
        title(titles{k}); grid on;
    end
end