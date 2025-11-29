clc
clear
% coupled penduli
% Marcello Farina, 19/12/2019, updated on 25/10/2020
seed = 123;
rng(seed)

k=0.02; %N/m
% k=2; %N/m
% k=200; %N/m

l = 1; %m
m = 1; %kg
g = 9.8; %m/s2
a = l;
h = 0.5; %s
N = 2; % Number of subsystems

A=[0 1 0 0
    g/l-k*a^2/m/l^2 0 k*a^2/m/l^2 0
    0 0 0 1
    k*a^2/m/l^2 0 g/l-k*a^2/m/l^2 0];
B=[[0 1/m/l^2 0 0]',[0 0 0 1/m/l^2]'];

C=eye(4);

%  1. Modelling
% ---------------------------------

% Decompose the state and input vectors into subvectors, consistently with the 
% physical description of the system. Obtain the corresponding decomposed model.

B1 = B(:,1);
B2 = B(:,2);
C1 = C(1:2,:);
C2 = C(3:end,:);

% Generate the system matrices (both continuous-time and discrete-time, the 
% latter with a sampling time selected compatibly with the continuous-time dynamics

% systemCT = ss(A,B,C,D);
% c2d(systemCT,h) does the same thing, and yields same results. 

% F = e^Ah
F = expm(A*h);
% G formula from solving the integral and applying the limits.
G = A \ (expm(A*h) - expm(A*0)) * B;
H = C;

G1 = G(:,1);
G2 = G(:,2);
H1 = H(1:2,:);
H2 = H(3:end,:);

% 2. Analysis - DT Centralized
% ---------------------------------

isStable(F);
ContStruc = ones(N,N);
discrete_centralized_fm = di_fixed_modes(F,{G1,G2},{H1,H2},2,ContStruc,4);
n_fm = size(discrete_centralized_fm);
n_fm = n_fm(1);
if (n_fm == 0)
    disp("The discrete centralized system has no fixed modes")
else
    disp("Fixed modes found at:")
    for i = 1:n_fm
        disp(discrete_centralized_fm(i))
    end
end

% 3. Design - LMI

% Generation of matrices 𝐿 and 𝑃 compatible with the information exchange structure
% sdpvar is used to define YALMIPs symbolic decision variables.
n_states = 4;
m_inputs = 2;
P=sdpvar(n_states);
L=sdpvar(m_inputs,n_states);

% Define the constraint LMI, solver can't handle zero so we use 1e-2
%  A. Stability - will provide a solution set (P,L), maybe not the best one, but definitely a stable one
LMIconstr_stability=[[P-F*P*F'-F*L'*G'-G*L*F' , G*L;
                                L'*G'         , P  ] >= 1e-2*eye(n_states*2)];

% Solve and simulate the simple stabitlity LMI

Kx = solveLMI(LMIconstr_stability,P,L);
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable.")
isStable(F_cl);
generate_plots(F_cl,h,n_states)

% B. to make all the eigenvalues of 𝐴 + 𝐵𝐾𝑥 satisfy |𝜆| < 𝜌.
%  This is possible by simply requiring that 1/𝜌 * (𝐴 +𝐵𝐾𝑥) is schur stable
%  the eigenvalues of 1/𝜌 * (𝐴 +𝐵𝐾𝑥) are 𝜆𝑖/𝜌

rho = 0.2; % |𝜆| < 0.2
LMIconstr_speed=[[rho^2*P-F*P*F'-F*L'*G'-G*L*F' , G*L;
                                L'*G'           , P  ] >= 1e-2*eye(n_states*2)];

% Solve and simulate the speedy LMI
Kx = solveLMI(LMIconstr_speed,P,L);
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable.")
isStable(F_cl);
generate_plots(F_cl,h,n_states)


% --- H2 Design: Define Weights ---
% J = \sum (x^T Q x + u^T R u)
Q = diag([1000,10,1000,10]); % High penaly on position error, low on speed error
R = eye(2); % Standard control effort penalty
% Noise on all 4 states
Gw = eye(4); 
%  4 states, 2 inputs --> z_k (6x1) =  H*x(4x1) + D_u*u(2x1)
%  H must be 6x4, D_u must be 6x2, so now we know how many zeros to pad
H = [
    sqrtm(Q);
    zeros(2,4)
    ];

Du = [
        zeros(4,2);
        sqrtm(R)
      ];


S = sdpvar(6, 6, 'symmetric'); % 6x6, size of z_k output (4 states + 2 inputs)

% --- H2 Design: Define LMI Constraints ---

% LMI 1: The Stability Constraint 
LMI_1 = [ [P - F*P*F' - F*L'*G' - G*L*F' - Gw*Gw',  G*L];
          [L'*G',                                   P] ] >= 1e-6 * eye(n_states*2);

% LMI 2: The Performance Constraint 
LMI_2 = [ [S,                           H*P + Du*L];
          [(H*P + Du*L)',     P] ] >= 1e-6 * eye(2*n_states+m_inputs);

% Combine all constraints
LMIconstr_H2 = [LMI_1, LMI_2];
options = sdpsettings('verbose', 0); % Define solver settings to be quiet
J = optimize(LMIconstr_H2,trace(S),options);
if J.problem 
    disp("Unfeasible")
end

L = double(L);
P = double(P);
Kx = L/P;
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable.")
isStable(F_cl);
generate_plots(F_cl,h,n_states)

% H-inf
gamma = sdpvar(1);
% --- H-Infinity Design Step 3: Define LMI Constraint ---
% We re-use H_perf and Du_perf
H_perf = H;
Du_perf = Du;

% We assume Dw = 0 (no direct noise-to-performance feedthrough)
Dw = zeros(size(H_perf,1), size(Gw,2)); % This will be 6x4

% Build the 4x4 Block LMI from the slides (Page 86)
% Block (1,1)
LMI_11 = P;
% Block (1,2)
LMI_12 = F*P + G*L;
% Block (1,3)
LMI_13 = Gw;
% Block (1,4)
LMI_14 = zeros(n_states, size(H_perf,1));

% Block (2,2)
LMI_22 = P;
% Block (2,3)
LMI_23 = zeros(n_states, size(Gw,2));
% Block (2,4)
LMI_24 = (H_perf*P + Du_perf*L)'; % Note the transpose

% Block (3,3)
LMI_33 = gamma * eye(size(Gw,2));
% Block (3,4)
LMI_34 = Dw'; % Note the transpose

% Block (4,4)
LMI_44 = gamma * eye(size(H_perf,1));

% Assemble the full matrix
LMI_Hinf = [ LMI_11, LMI_12, LMI_13, LMI_14;
             LMI_12',LMI_22, LMI_23, LMI_24;
             LMI_13',LMI_23',LMI_33, LMI_34;
             LMI_14',LMI_24',LMI_34',LMI_44 ];

constraints = [LMI_Hinf >= 1e-6 * eye(size(LMI_Hinf,1)), P >= 1e-6 * eye(n_states)];

% --- H-Infinity Design Step 4: Solve and Extract ---
disp('Solving H-Infinity LMI... (This may take a moment)');
% We tell it to MINIMIZE gamma
sol = optimize(constraints, gamma, options);

if sol.problem == 0
    % Success!
    K_hinf = double(L) / double(P);
    
    disp('H-Infinity Controller (K_hinf) successfully designed:');
    disp(K_hinf);
    
    fprintf('Optimal H-Infinity Gain (gamma): %.4f\n', double(gamma));
    
    % Verify and Plot
    F_cl_hinf = F + G * K_hinf;
    rho_hinf = max(abs(eig(F_cl_hinf)));
    
    fprintf('Closed-loop spectral radius: %.4f\n', rho_hinf);
    generate_plots(F_cl_hinf, h, n_states);
else
    % Failure
    disp('!!! H-Infinity LMI was Infeasible !!!');
    disp(sol.info);
end


function Kx = solveLMI(LMIconstr,P,L)
    options = sdpsettings('verbose', 0); % Define solver settings to be quiet
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
    % Compute the eigenvalues and the spectral radius of the (discrete-time) system. 
    % Is it open-loop asymptotically stable?
    rho = max(abs(eig(F)));
    if (rho > 1) 
        fprintf ('the absolute spectral radius of the discrete system is %.2f > 1 and is therfore unstable.\n',rho)
    else
        fprintf("The system is stable, spectral radius = |%.2f| < 1\n", rho)
    end
end

function generate_plots(F_cl,h,n_states)
    % 1. DEFINE SIMULATION PARAMETERS
    
    % Pick a non-zero initial condition
    x0_1 = (pi/6); % Initial angle for pendulum 1 (3 deg)
    x0_2 = 0;
    x0_3 = -(pi/6); % Initial angle for pendulum 2 (3 deg)
    x0_4 = 0;
    
    x0 = [x0_1; x0_2; x0_3; x0_4];   
    
    n_steps = 10; % Number of steps to simulate
    t_final = h * n_steps;
    t = 0:h:t_final; % Create a time vector for plotting
    
    % n_states is from the LMI code 
    x_history = zeros(n_states, n_steps + 1);
    x_history(:, 1) = x0; 
    
    % 4. RUN THE SIMULATION LOOP
    current_x = x0;
    for k = 1:n_steps
        % This is the core simulation step
        next_x = F_cl * current_x;
        
        % Store the result
        x_history(:, k+1) = next_x;
        
        % Update for the next loop
        current_x = next_x;
    end

    figure;
    sgtitle('Closed-Loop State Trajectories'); 
    
    % Subplot 1: x1 (Position 1)
    subplot(2, 2, 1);
    plot(t, x_history(1, :));
    title('State x1 (Position 1)');
    xlabel('Time (s)');
    ylabel('Position');
    grid on;
    
    % Subplot 2: x2 (Velocity 1)
    subplot(2, 2, 2);
    plot(t, x_history(2, :));
    title('State x2 (Velocity 1)');
    xlabel('Time (s)');
    ylabel('Velocity');
    grid on;
    
    % Subplot 3: x3 (Position 2)
    subplot(2, 2, 3);
    plot(t, x_history(3, :));
    title('State x3 (Position 2)');
    xlabel('Time (s)');
    ylabel('Position');
    grid on;
    
    % Subplot 4: x4 (Velocity 2)
    subplot(2, 2, 4);
    plot(t, x_history(4, :));
    title('State x4 (Velocity 2)');
    xlabel('Time (s)');
    ylabel('Velocity');
    grid on;
end


