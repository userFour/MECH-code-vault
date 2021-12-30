% -----------------------------------------------------------------------------
% MECH463 Car Project
%
%
% Author: Markos    Date: 2021-11-23
% -----------------------------------------------------------------------------
%
% -----------------------------------------------------------------------------
%
% Sample calls: Just run it lol
% R = @(t) (0.001*sin(0.05*t) + 0.005*sin(0.5*t) + ...
%     0.001*sin(10*t) + 0.002*sin(50*t) + 0.0001*sin(100*t) + ...
%     0.0003*sin(250*t) +  0.004*sin(500*t));
%
% plot(t, R(t), 'r-')
%
% x_init = [0.05; -0.05; -0.05; 0; 0; 0; 0];
% v_init = [-0.1; 0.1; 0.1; 0; 0; 0; 0];
%
% x_init = [0; 0; 0; 0; 0; 0; 0];
% v_init = [0; 0; 0; 0; 0; 0; 0];
%
% -----------------------------------------------------------------------------
%

clc; clear; close;
fprintf("Start program!");

% Physical Paraments
m_car = 1000;   % [kg]
l_car = 4.7;    % [m]
w_car = 1.7;    % [m]
h_car = 1.7;    % [m]
c = 1000;       % [Ns/m]
kf = 10000;     % [N/m]
kr = 10000;     % [N/m]
kt = 10000;     % [N/m]
m_t = 12;       % [kg]

% Computed physical parameters
a1 = 0.6*l_car;
a2 = l_car - a1;
b1 = 0.6*w_car;
b2 = w_car - b1;
I_y = (1/12) * m_car * (w_car^2 + h_car^2);
I_z = (1/12) * m_car * (l_car^2 + h_car^2);
fprintf('\nI_y = %f', I_y);
fprintf('\nI_z = %f', I_z);

% Initial conditions
% x, phi, theta, xa, xb, xc, xd
x_init = [0; 0; 0; 0; 0; 0; 0];
v_init = [0; 0; 0; 0; 0; 0; 0];

% Mass, damping, stiffness matrix
M = [m_car, 0, 0, 0, 0, 0, 0; ...
     0, I_y, 0, 0, 0, 0, 0; ...
     0, 0, I_z, 0, 0, 0, 0; ...
     0, 0, 0, m_t, 0, 0, 0; ...
     0, 0, 0, 0, m_t, 0, 0; ...
     0, 0, 0, 0, 0, m_t, 0; ...
     0, 0, 0, 0, 0, 0, m_t];

C = [4*c, 2*c*(a2 - a1), 2*c*(b2 - b1), -c, -c, -c, -c; ...
     2*c*(a2 - a1), 2*c*(a2^2 + a1^2), c*(b1*a2 + b2*a1 - b1*a1 - b2*a2), c*a1, c*a1, -c*a2, -c*a2; ...
     2*c*(b1 - b2), c*(a2*b1 + a1*b2 - a1*b1 - a2*b2), 2*c*(b1^2 + b2^2), c*b2, -c*b1, -c*b1, c*b2; ...
     -c, c*a1, c*b2, c, 0, 0, 0; ...
     -c, c*a1, -c*b1, 0, c, 0, 0; ...
     -c, -c*a2, -c*b1, 0, 0, c, 0; ...
     -c, -c*a2, c*b2, 0, 0, 0, c];

K = [2*(kf + kr), 2*(kr*a2 - kf*a1), (kf + kr)*(b1 - b2), -kf, -kf, -kr, -kr; ...
     2*(kr*a2 - kf*a1), 2*(kf*a1^2 + kr*a2^2), (kr*b1*a2 + kf*b2*a1 - kr*b2*a2 - kf*b1*a1), kf*a1, kf*a1, -kr*a2, -kr*a2; ...
     (kf + kr)*(b1 - b2), (kr*a2*b1 + kf*a1*b1 - kf*a1*b1 - kr*a2*b2), (kf + kr)*(b1^2 + b2^2), kf*b2, -kf*b1, -kr*b1, kr*b2; ...
     -kf, kf*a1, kf*b2, (kt + kf), 0, 0, 0; ...
     -kf, kf*a1, -kf*b1, 0, (kt + kf), 0, 0; ...
     -kr, -kr*a2, -kr*b1, 0, 0, (kt + kr), 0; ...
     -kr, -kr*a2, kr*b2, 0, 0, 0, (kt + kr)];

% Cast the system into state space
A_ss = [zeros(7, 7), eye(7, 7); -inv(M)*K, -inv(M)*C];
B_ss = [zeros(7, 7); inv(M)];
C_ss = [eye(7, 7), zeros(7, 7)];
D_ss = zeros(7, 7);

% Generate the state space model
ss_model = ss(A_ss, B_ss, C_ss, D_ss);

% Define time range, intitial conditions
t = linspace(0, 20, 10E3);
x_0 = [x_init; v_init];

% Set up the forcing function
base_motion = (0.001*sin(0.05*t) + 0.005*sin(0.5*t) + ...
          0.001*sin(10*t) + 0.002*sin(50*t) + 0.0001*sin(100*t) + ...
          0.0003*sin(250*t) +  0.004*sin(500*t));

forcing_fxn = kt * base_motion;

F_t = [zeros(1,length(t)); ...
       zeros(1,length(t)); ...
       zeros(1,length(t)); ...
       forcing_fxn; ...
       forcing_fxn; ...
       forcing_fxn; ...
       forcing_fxn];

% Solve for output
x = lsim(ss_model, F_t, t, x_0);

% Plot results
hold on;
subplot(3, 2, 1);
plot(t, x(:, 1), 'b', 'LineWidth', 1.5); % Chassis
%plot(t, base_motion, 'k'); % Plot the input
xlabel("Time [s]");
ylabel("Position [m]");
title("Chassis Vibrations");

subplot(3, 2, 2);
plot(t, x(:, 2), 'r', 'LineWidth', 1.5); % Phi
%plot(t, base_motion, 'k'); % Plot the input
xlabel("Time [s]");
ylabel("Position [m]");
title("Phi");

subplot(3, 2, 3);
plot(t, x(:, 3), 'g', 'LineWidth', 1.5); % Theta
%plot(t, base_motion, 'k'); % Plot the input
xlabel("Time [s]");
ylabel("Position [m]");
title("Theta");

subplot(3, 2, 4);
plot(t, x(:, 4), 'c', 'LineWidth', 1.5); % Wheel A
plot(t, x(:, 5), 'c', 'LineWidth', 1.5); % Wheel B
plot(t, x(:, 6), 'c', 'LineWidth', 1.5); % Wheel C
plot(t, x(:, 7), 'c', 'LineWidth', 1.5); % Wheel D
%plot(t, base_motion, 'k'); % Plot the input
xlabel("Time [s]");
ylabel("Position [m]");
title("Wheels");

subplot(3, 2, [5, 6]);
plot(t, base_motion, 'k', 'LineWidth', 1.5); % Road
xlabel("Time [s]");
ylabel("Position [m]");
title("Road");

grid on;
hold off;

% Determine max magnitude of chassis displacement
fprintf("\nMax displacement: %0.4f [mm]", (max(x(:, 1))*1000))

fprintf("\nEnd Program!");

