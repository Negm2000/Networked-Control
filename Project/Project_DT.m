clc
clear
close all

% coupled penduli
% Marcello Farina, 19/12/2019, updated on 25/10/2020

k = 0.2;   % N/m
% k = 2;   % N/m
% k = 200; % N/m
l = 1;   % m
m = 1;   % kg
g = 9.8; % m/s2
a = l;
h = 0.5; % s
N = 2;   % Number of subsystems
A = [
    0               1       0         0
    g/l-k*a^2/m/l^2 0 k*a^2/m/l^2     0
    0               0      0          1
    k*a^2/m/l^2     0 g/l-k*a^2/m/l^2 0
    ];
B=[[0 1/m/l^2 0 0]',[0 0 0 1/m/l^2]'];
C=eye(4);

%  1. Modelling
% ---------------------------------
% Decompose the state and input vectors into subvectors
B1 = B(:,1);
B2 = B(:,2);
C1 = C(1:2,:);
C2 = C(3:end,:);

% Generate the system matrices 
F = expm(A*h);
% G formula 
G = A \ (expm(A*h) - expm(A*0)) * B;
H = C;
G1 = G(:,1);
G2 = G(:,2);
H1 = H(1:2,:);
H2 = H(3:end,:);

% 2. Analysis - DT Centralized
% ---------------------------------
isStable(F);

% 3. Design - LMI
n_states = 4;
m_inputs = 2;
P=sdpvar(n_states);
L=sdpvar(m_inputs,n_states);

%  A. Stability
% ---------------------------------
LMIconstr_stability=[[P-F*P*F'-F*L'*G'-G*L*F' , G*L;
                                L'*G'         , P  ] >= 1e-2*eye(n_states*2)];
Kx = solveLMI(LMIconstr_stability,P,L);
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable (Stability LMI).")
isStable(F_cl);

% VISUALIZATION A
generate_plots_smooth(Kx, A, B, h, n_states, 'Stability LMI');
animate_pendulums(Kx, A, B, h, 'Stability Animation'); 

%  B. Pole Placement (Speed)
% ---------------------------------
rho = 0.12; % |𝜆| < 0.2
LMIconstr_speed=[[rho^2*P-F*P*F'-F*L'*G'-G*L*F' , G*L;
                                L'*G'           , P  ] >= 1e-2*eye(n_states*2)];
Kx = solveLMI(LMIconstr_speed,P,L);
F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable (Speed LMI).")
isStable(F_cl);

% VISUALIZATION B
generate_plots_smooth(Kx, A, B, h, n_states, 'Speed LMI');
animate_pendulums(Kx, A, B, h, 'Speed LMI Animation'); 

%  C. H2 Design
% ---------------------------------
P=sdpvar(n_states); 
L=sdpvar(m_inputs,n_states);

% Bryson's Rule
%Q_{ii} = 1/(max acceptable value for x_i)^2
%R_{ii} = 1/(max acceptable value for u_i)^2
Q = diag([0.01^-2, 0.1^-2, 0.01^-2, 0.1^-2]); 
R = 10*eye(2); 
Gw = eye(4); 
H = [sqrtm(Q); zeros(2,4)];
Du = [zeros(4,2); sqrtm(R)];
S = sdpvar(6, 6, 'symmetric'); 

LMI_1 = [ [P - F*P*F' - F*L'*G' - G*L*F' - Gw*Gw',  G*L];
          [L'*G',                                   P] ] >= 1e-6 * eye(n_states*2);

LMI_2 = [ [S,                           H*P + Du*L];
          [(H*P + Du*L)',     P] ] >= 1e-6 * eye(2*n_states+m_inputs);

LMIconstr_H2 = [LMI_1, LMI_2];
options = sdpsettings('verbose', 0); 
J = optimize(LMIconstr_H2,trace(S),options);
if J.problem 
    disp("Unfeasible")
end

% Manual extraction (as in your code)
L_val = double(L);
P_val = double(P);
Kx = L_val/P_val;

F_cl = F + G*Kx;
disp("Check if the centralized discrete system is stable (H2).")
isStable(F_cl);

% VISUALIZATION C
generate_plots_smooth(Kx, A, B, h, n_states, 'H2 Design');
animate_pendulums(Kx, A, B, h, 'H2 Animation'); 

%  D. H-Infinity Design
% ---------------------------------
P = sdpvar(n_states);
L = sdpvar(m_inputs,n_states);

gamma = sdpvar(1);
H_perf = H;
Du_perf = Du;
Dw = zeros(size(H_perf,1), size(Gw,2)); 

% Build the 4x4 Block LMI 
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
sol = optimize(constraints, gamma, options);

if sol.problem == 0
    K_hinf = double(L) / double(P);
    disp('H-Infinity Controller (K_hinf) successfully designed.');
    fprintf('Optimal H-Infinity Gain (gamma): %.4f\n', double(gamma));
    
    F_cl_hinf = F + G * K_hinf;
    rho_hinf = max(abs(eig(F_cl_hinf)));
    fprintf('Closed-loop spectral radius: %.4f\n', rho_hinf);
    
    % VISUALIZATION D
    generate_plots_smooth(K_hinf, A, B, h, n_states, 'H-Inf Design');
    animate_pendulums(K_hinf, A, B, h, 'H-Inf Animation');
else
    disp('H-Infinity LMI was Infeasible');
    disp(sol.info);
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
    rho = max(abs(eig(F)));
    if (rho > 1) 
        fprintf ('The absolute spectral radius is %.2f > 1 (Unstable).\n',rho)
    else
        fprintf("The system is stable, spectral radius = |%.2f| < 1\n", rho)
    end
end

% --- HYBRID SMOOTH PLOTTING ---
function generate_plots_smooth(K, A, B, h, n_states, plotTitle)
    % Physics simulation settings
    dt_physics = 0.01; 
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
        
        % Discrete Control Update
        if t_now >= next_sample_time
            if ~isempty(K)
                u_current = K * current_x; 
            else
                u_current = [0;0];
            end
            next_sample_time = next_sample_time + h;
        end
        
        % Continuous Physics Update
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

% --- ANIMATION ---
function animate_pendulums(K, A, B, h, plotTitle)
    dt_physics = 0.05; 
    n_steps = 15; % Longer duration for animation
    t_final = h * n_steps;
    t = 0:dt_physics:t_final;
    
    x0 = [pi/6; 0; -pi/6; 0];
    current_x = x0;
    
    u_current = [0;0];
    next_sample_time = 0;
    
    figure('Name', plotTitle, 'Color', 'w');
    l = 1; 
    offset = 2; 
    
    plot([-1, offset+1], [0,0], 'k-', 'LineWidth', 2); hold on;
    
    rod1 = plot([0,0], [0,0], 'b-', 'LineWidth', 2);
    bob1 = plot(0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    rod2 = plot([0,0], [0,0], 'r-', 'LineWidth', 2);
    bob2 = plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    spring = plot([0,0], [0,0], 'g--', 'LineWidth', 1);
    
    axis([-1.5, offset+1.5, -1.5, 1.5]); % Adjusted Y axis for inverted potential
    axis equal; grid on;
    title(['Animation: ' plotTitle]);
    
    disp('Starting animation...');
    pause(0.5);
    
    for i = 1:length(t)
        if t(i) >= next_sample_time
            if ~isempty(K), u_current = K * current_x; end
            next_sample_time = next_sample_time + h;
        end
        
        dx = A * current_x + B * u_current;
        current_x = current_x + dx * dt_physics;
        
        theta1 = current_x(1);
        theta2 = current_x(3);
        
        % Inverted logic (matches your positive g/l)
        x1 = l * sin(theta1);
        y1 = l * cos(theta1); 
        
        x2 = offset + l * sin(theta2);
        y2 = l * cos(theta2);
        
        set(rod1, 'XData', [0, x1], 'YData', [0, y1]);
        set(bob1, 'XData', x1, 'YData', y1);
        set(rod2, 'XData', [offset, x2], 'YData', [0, y2]);
        set(bob2, 'XData', x2, 'YData', y2);
        set(spring, 'XData', [x1, x2], 'YData', [y1, y2]);
        
        drawnow; 
        pause(0.02);
    end
end