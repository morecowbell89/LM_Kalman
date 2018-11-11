% LM KF Class
% HW2 Problem 2
% 2018-11-10
% Implement Kalman Filter for Gyro

%% 
clear all;
close all;

% sim time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.1;
t_start = 0;
t_end = 1000;
time = t_start:dt:t_end;

% truth parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_nu = 3e-7; % Angle random walk (rad/sec^0.5)
sigma_u = 3e-10; % rate random walk (rad/sec^(3/2))
sigma_S = 17e-6; % noise on measured angle (rad)
omega_in = 1e-7; % constant rate truth (rad/sec)
Phi = [1 0; 0 1]; % calculate state transition model; states are [theta beta]
Lambda = [dt 0]'; % used for inputing omega*dt into theta
Gamma = [0 sigma_u*sqrt(dt)]'; 
H = [1 0];

% kalman filter parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = [(sigma_nu^2*dt + (1/3)*sigma_u^2*dt^3) -0.5*sigma_u^2*dt;
    -0.5*sigma_u^2*dt sigma_u^2*dt]; % From slides, process noise covari. 
Gamma_kf = [dt 0; 0 0];
R = sigma_S^2; % Measurment noise covariance
Phi_kf = [1 -dt; 0 1]; 
H_kf = H; 

% preallocating vectors
%%%%%%%%%%%%%%%%%%%%%%%%%
omega = zeros(1,length(time)+1);
P_prior = zeros(2,2,length(time)+1);
x = zeros(2,length(time)+1);
x_hat_prior = zeros(2,length(time)+1);
y = zeros(1,length(time));
e_meas = zeros(1,length(time));
x_hat_post = zeros(2,length(time));
P_post = zeros(2,2,length(time));
e_states = zeros(2,length(time));
one_sigma_theta = zeros(1,length(time));
one_sigma_beta = zeros(1,length(time));
K = zeros(2,length(time));

%% Propogate
for k = 1:length(time)
    if k == 1 
        % Initialize states
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initalize state estimate
        x_hat_prior(:,k) = [0 0]'; 
        % initialize covariance estimate
        P_prior(:,:,k) = (1e-4)*eye(2); % P_initial
        % initialize truth state
        x(:,k) = [0 1e-6]'; % truth initial state [theta beta]'
    end
    
    % Propogate truth states
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % propogate true input
    % omega_in constant, setting both omega(k) and omega(k+1) to omega_in
    omega(k) = omega_in; 
    omega(k+1) = omega_in; 
    
    % propogate truth states Theta and Beta, x(1) = theta, x(2) = beta
    x(:,k+1) = Phi*x(:,k) + Lambda*omega(k+1) + Gamma*randn;
    % calculate K gain
    K(:,k) = P_prior(:,:,k)*H_kf'*inv(H_kf*P_prior(:,:,k)*H_kf' + R);
    
    % Take in measurments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute theta measured
    y(:,k) = H*x(:,k) + sigma_S*randn;
    % compute omega measured
    omega_meas = omega(k+1) +...
        0.5*(x(2,k+1)+x(2,k)) +...
        randn*sqrt(sigma_nu^2/dt + (sigma_u^2*dt)/12);
    
    % After measurment steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate residual from measurment and estimate of state prior to
    % measurement
    e_meas(:,k) = y(:,k) - H_kf*x_hat_prior(:,k);    
    % update state posteriori
    x_hat_post(:,k) = x_hat_prior(:,k) + K(:,k)*e_meas(:,k);
    % update covariance posteriori
    P_post(:,:,k) = (eye(2) - K(:,k)*H_kf)*P_prior(:,:,k);    
    e_states(:,k) = x(:,k) - x_hat_post(:,k); % record for plotting    
    % Propogate estimated state and covariance forward time step
    % propogate estimates forward in time
    x_hat_prior(:,k+1) = Phi_kf*x_hat_post(:,k) + Lambda*omega_meas;
    % Propogate estimated covariance forward time step
    P_prior(:,:,k+1) = Phi_kf*P_post(:,:,k)*Phi_kf' + Gamma_kf*Q*Gamma_kf';
    
    % calculate one_sigma values for theta and beta from covariance matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    one_sigma_theta(k) = sqrt(P_post(1,1,k));
    one_sigma_beta(k) = sqrt(P_post(2,2,k));
    
    
end

%% Plots

% Theta and Beta Gains
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
ax1 = subplot(2,1,1);
plot(time,K(1,:));
ylabel('Theta Gain')
ax2 = subplot(2,1,2);
plot(time,K(2,:));
ylabel('Beta Gain')
linkaxes([ax1,ax2],'x')

% Theta and Beta Errors with 3-sig values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
% theta errors and 3-sigma
ax1 = subplot(2,1,1);
plot(time,e_states(1,:));
hold on;
plot(time,3*one_sigma_theta,'r');
plot(time,-3*one_sigma_theta,'r');
xlabel('Time (s)');
ylabel('\theta Error (rad)');
legend({'\theta','\theta 3-\sigma'});

% beta errors and 3-sigma
ax2 = subplot(2,1,2);
plot(time,e_states(2,:));
hold on;
plot(time,3*one_sigma_beta,'r');
plot(time,-3*one_sigma_beta,'r');
legend('\beta','\beta 3-\sigma');
xlabel('Time (s)');
ylabel('\beta Error (rad/sec)')
linkaxes([ax1,ax2],'x')

% Theta and Beta Estimates with Theta and Beta True
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()

% theta estimate and theta true
ax1 = subplot(2,1,1);
plot(time,x(1,1:(end-1)));
hold on;
plot(time,x_hat_post(1,:),'r');
ylabel('\theta (rad)');
legend({'\theta True','\theta Estimate'});
% beta estimate and beta true
ax2 = subplot(2,1,2);
plot(time,x(2,1:(end-1)));
hold on;
plot(time,x_hat_post(2,:),'r');
legend({'\beta True','\beta Estimate'});
ylabel('\beta (rad/s)')
linkaxes([ax1,ax2],'x')

%% Bonus: Compare Sim Covariance to Analytical
P_the_bet = -sigma_u*sigma_S*sqrt(dt);
P_the_the = sigma_S*sqrt((sigma_nu^2-2*P_the_bet)*dt);
P_bet_bet = sigma_u*sqrt(sigma_nu^2-2*P_the_bet);
P_an = [P_the_the P_the_bet;
        P_the_bet P_bet_bet];
P_est_ss = P_post(:,:,end);

disp(P_an - P_est_ss); % Results agree to ~1e-11
%{
   1.0e-11 *

    0.1436   -0.0001
   -0.0001    0.0000
%}

