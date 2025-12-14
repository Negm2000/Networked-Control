clc
clear
close all

% Project 7: Coupled Penduli

k = 200;   % N/m (Spring stiffness)
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
% Unstable if any eigenvalue magnitude > 1
disp("Check if the centralized discrete system is stable (Open Loop).")
isStable(F);
disp("---------------------------------");
disp("Controller Design With LMIs");
disp("---------------------------------");
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
analyze_system(F_cl, h, Kx, n_states, 'Stability LMI');

%  B. Pole Placement (Spectral Radius / Speed)
% ---------------------------------
% Goal: Force eigenvalues |lambda| < rho
% How: Require (1/rho)*(F+GK) to be Schur stable.
rho = 0.61;
LMIconstr_speed=[[rho^2*P-F*P*F'-F*L'*G'-G*L*F' , G*L;
    L'*G'           , P  ] >= 1e-2*eye(n_states*2)];

Kx = solveLMI(LMIconstr_speed,P,L);
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable (Speed LMI).")
isStable(F_cl);
% VISUALIZATION B
analyze_system(F_cl, h, Kx, n_states, 'Speed LMI');

%  C. H2 Optimal Design
% ---------------------------------
% Goal: Minimize Trace(S), where S is an upper bound on output energy z_k.
P=sdpvar(n_states);
L=sdpvar(m_inputs,n_states);

% Bryson's Rule for Weights
% Q_{ii} = (max acceptable value for x_i)^-2
% R_{ii} = (max acceptable value for u_i)^-2
Q = 100*eye(4);
R = 0.1*eye(2); % Input penalty

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
analyze_system(F_cl, h, Kx, n_states, 'H2 Design');

%  D. H-Infinity Robust Design
% ---------------------------------
% Goal: Minimize worst-case gain gamma from disturbance w to output z.
% Uses the Bounded Real Lemma in LMI form.

P = sdpvar(n_states);
L = sdpvar(m_inputs,n_states);
gamma = sdpvar(1);

Dw = zeros(size(H,1), size(Gw,2));

% Build the 4x4 Block Matrix for Bounded Real Lemma (Page 86)
LMI_11 = P;
LMI_12 = F*P + G*L;
LMI_13 = Gw;
LMI_14 = zeros(n_states, size(H,1));

LMI_22 = P;
LMI_23 = zeros(n_states, size(Gw,2));
LMI_24 = (H*P + Du*L)';

LMI_33 = gamma * eye(size(Gw,2));
LMI_34 = Dw';

LMI_44 = gamma * eye(size(H,1));

LMI_Hinf = [ LMI_11, LMI_12, LMI_13, LMI_14;
    LMI_12',LMI_22, LMI_23, LMI_24;
    LMI_13',LMI_23',LMI_33, LMI_34;
    LMI_14',LMI_24',LMI_34',LMI_44 ];

constraints = [LMI_Hinf >= 1e-6 * eye(size(LMI_Hinf,1)), P >= 1e-6 * eye(n_states)];

sol = optimize(constraints, gamma, options); % Minimize gamma

if sol.problem == 0
    K_hinf = double(L) / double(P);
    F_cl_hinf = F + G * K_hinf;
    disp("Check if the centralized discrete system is stable (H-inf).")
    isStable(F_cl_hinf);
    fprintf('Optimal H-Infinity Gain (gamma): %.4f\n', double(gamma));
    % VISUALIZATION D
    analyze_system(F_cl_hinf, h, K_hinf, n_states, 'H-Inf Design');
else
    disp('H-Infinity LMI was Infeasible');
    disp(sol.info);
end

%  E. Damping Design
% -----------------------------------------------------------------
% Goal: Force eigenvalues to be Positive Real to prevent ringing.
% Method: Shifted Circle LMI. Center (alpha) = 0.4, Radius (r) = 0.4.
P = sdpvar(n_states);
L = sdpvar(m_inputs, n_states);

alpha = 0.4;   % Center of the circle
r_damp = 0.4;  % Radius of the circle

% We simply shift the A matrix (F) by alpha*Identity
F_shifted = F - alpha * eye(n_states);

% LMI: Schur Stability on Shifted Matrix
% Condition: (F_shift + G*K) * P * (F_shift + G*K)' < r^2 * P
LMIconstr_damping = [[r_damp^2*P - F_shifted*P*F_shifted' - F_shifted*L'*G' - G*L*F_shifted', G*L;
    L'*G'                                                      , P ] >= 1e-6*eye(n_states*2)];

Kx = solveLMI(LMIconstr_damping,P,L);
F_cl = F + G*Kx;

disp("Check if the centralized discrete system is stable (Damping LMI).")
isStable(F_cl);

% VISUALIZATION E
analyze_system(F_cl, h, Kx, n_states, 'Damping Design');

%  F. Multi-Objective Design
% -----------------------------------------------------------------
% Combines:
% 1. Speed (Pole Placement)
% 2. Damping (Cardioid)
% 3. Effort Minimization using variables Kl and Ky

disp("---------------------------------");
disp("Part F: Multi-Objective Design");
disp("---------------------------------");

% 1. Define Variables
P = sdpvar(n_states, n_states, 'symmetric');
L = sdpvar(m_inputs, n_states, 'full');
Kl = sdpvar(1, 1); % Scalar proxy for norm(L)
Ky = sdpvar(1, 1); % Scalar proxy for norm(inv(P))

% 2. Tuning Parameters
% Relaxed speed constraint to allow the solver "breathing room"
% for the control effort minimization.
rho_speed = 0.35;  % Decay rate (Close to 1 is easier, Close to 0 is faster)

% Damping (Cardioid) parameters
% We keep the region large to ensure feasibility
alpha = 0.3;
r_damp = 0.3;
F_shifted = F - alpha * eye(n_states);

% 3. Define Constraints

% Constraint A: Speed / Stability
% (rho^2 * P) - (F*P + G*L) * P^-1 * (F*P + G*L)' > 0
LMI_Speed = [[rho_speed^2*P,       (F*P + G*L)';
    (F*P + G*L),         P           ] >= 1e-6*eye(2*n_states)];

% Constraint B: Damping (Shifted Circle)
% Forces poles away from the unit circle edge (reduces ringing)
LMI_Damping = [[r_damp^2*P,        (F_shifted*P + G*L)';
    (F_shifted*P + G*L), P                  ] >= 1e-6*eye(2*n_states)];

% Constraint C: Control Effort ("Kl" Logic)
% [Kl*I, L'; L, I] > 0  ==>  Kl > L*L' (Bounds norm of L)
LMI_Kl = [[Kl*eye(n_states), L';
    L,                eye(m_inputs)] >= 1e-6*eye(n_states + m_inputs)];

% Constraint D: Conditioning ("Ky" Logic)
% [Ky*I, I; I, P] > 0   ==>  Ky*P > I  ==> P > 1/Ky * I (Bounds P away from 0)
% This prevents the solver from making P near-singular to cheat the speed constraint.
LMI_Ky = [[Ky*eye(n_states), eye(n_states);
    eye(n_states),    P            ] >= 1e-6*eye(2*n_states)];

% Combine Constraints
Constraints = [LMI_Speed, LMI_Damping, LMI_Kl, LMI_Ky, P >= 1e-6*eye(n_states)];

% 4. Optimization
% Cost Function: Weighted sum of Controller Gain (Kl) and Conditioning (Ky)
% Weighting Ky heavily forces P to be "well conditioned" (numerically stable)
J = 1*Kl + 100*Ky;

options = sdpsettings('verbose', 0, 'solver', 'sdpt3');

sol = optimize(Constraints, J, options);

if sol.problem == 0
    L_val = double(L);
    P_val = double(P);
    K_multi = L_val / P_val;

    disp('Feasible Solution Found!');
    fprintf('Minimization Result: Kl=%.2f, Ky=%.2f\n', double(Kl), double(Ky));

    % Verification
    F_cl_multi = F + G * K_multi;
    disp("Check stability (Multi-Objective).")
    isStable(F_cl_multi);

    % VISUALIZATION F
    analyze_system(F_cl_multi, h, K_multi, n_states, 'Multi-Objective');
else
    disp('LMI Infeasible.');
    disp(sol.info);
    disp('Try increasing rho_speed (closer to 1) or r_damp.');
end
% =========================================================================
% FUNCTIONS
% =========================================================================

function Kx = solveLMI(LMIconstr,P,L)
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
eigenvals = eig(F);
[rho, rho_idx] = max(abs(eigenvals));
if (rho >= 1)
    fprintf ('The spectral radius is |%.2f| >= 1 (Unstable).\n',eigenvals(rho_idx));
else
    fprintf("The spectral radius is  |%.2f| < 1 (Stable).\n", eigenvals(rho_idx));
end
end

function analyze_system(F_cl, h, K, n_states, plotTitle)
% 1. Simulation
[t, x_hist, u_hist] = simulate_discrete(F_cl, h, K, n_states);

% 2. Plotting
plot_simulation(t, x_hist, plotTitle);

% 3. Metrics
calculate_metrics(x_hist, u_hist, K, F_cl, h, t);
end

function [t, x_hist, u_hist] = simulate_discrete(F_cl, h, K, n_states)
n_steps = 50;
t = 0:h:(h*n_steps);

x0 = [pi/6; 0; -pi/6; 0];
x_hist = zeros(n_states, length(t));
x_hist(:, 1) = x0;

for k = 1:n_steps
    x_hist(:, k+1) = F_cl * x_hist(:, k);
end
u_hist = K * x_hist;
end

function plot_simulation(t, x_hist, plotTitle)
figure('Name', plotTitle);
sgtitle([plotTitle ' (Discrete Simulation)']);

titles = {'Pos 1 (x1)', 'Vel 1 (x2)', 'Pos 2 (x3)', 'Vel 2 (x4)'};
for k = 1:4
    subplot(2, 2, k);
    plot(t, x_hist(k, :), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    title(titles{k});
    grid on;
    xlabel('Time (s)');
end
end

function calculate_metrics(x_hist, u_hist, K, F_cl, h, t)
% Measured settling time 2% criterion
% equivalent to sqrt(theta1^2 + theta1_dot^2 + theta2^2 + theta2_dot^2)
traj_norm = vecnorm(x_hist);
threshold = 0.02 * max(traj_norm);
unsettled_indices = find(traj_norm > threshold);

if isempty(unsettled_indices)
    Ts = 0;
else
    Ts = t(unsettled_indices(end));
end

% Performance Metrics
% First sum sums over all inputs, second sum sums over all time
cont_energy = sum(sum(u_hist.^2)) * h;
state_energy = sum(sum(x_hist.^2)) * h;

% Theoretical Settling Time
rho = max(abs(eig(F_cl)));
if rho < 1
    Ts_theory = (log(0.02) / log(rho)) * h;
else
    Ts_theory = Inf;
end

fprintf("- Theoretical Max Settling Time (2%%) = %.2f\n", Ts_theory);
fprintf("- Measured Settling Time (2%%) (Ts) = %.2f\n", Ts);
fprintf("- Controller Gain Norm |K| = %.2f\n", norm(K));
fprintf("- Simulation Control Energy (Ju) = %.2f\n", cont_energy);
fprintf("- Simulation State Energy (Jx) = %.2f\n\n", state_energy);
end