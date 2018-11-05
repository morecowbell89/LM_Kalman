%% Problem 2
clear all;
close all;
clc;

% gyro error statistics
sigma_nu = 3e-7; %[rad/sec^0.5]
sigma_u = 3e-10; %[rad/sec^1.5]
sigma_S = 17e-6; % [rad]
% truth rate
theta_dot = 1e-7; %[rad/sec], spacecraft rate, constant
% sim timestep, also timestep of observer
dt = 0.1;
% observer gains
K_theta = 1;
K_beta = -1/dt;

% start simming
t_stop = 3000; % [sec]
t = 0:dt:t_stop;
% initialize
step_cnt = 1;

% truth initialization
beta(step_cnt) = 1e-6; % [rad/sec]
theta(step_cnt) = 0; % initial angle

% observer initialization
theta_hat_pri = 0; % [rad]
beta_hat_pri = 0; %[rad/sec]
% measurements
theta_meas(step_cnt) = theta(step_cnt) + sigma_S*randn();
omega_meas(step_cnt) = theta_dot + beta(step_cnt) + sigma_nu*randn();

% run observer from error measurments
e_meas = theta_meas(step_cnt) - theta_hat_pri;
theta_hat_post(step_cnt) = theta_hat_pri + K_theta*e_meas;
beta_hat_post(step_cnt) = beta_hat_pri + K_beta*e_meas;

e_truth(step_cnt) = theta_hat_post(step_cnt) - theta(step_cnt);

% start propogating
for step_cnt = 2:length(t)
    % propogate truth model first
    theta(step_cnt) = theta(step_cnt-1) + theta_dot*dt;
    beta(step_cnt) = beta(step_cnt-1) + sigma_nu*sqrt(dt)*randn();
    
    % gyro measure model
    theta_meas(step_cnt) = theta(step_cnt) + sigma_S*randn();
    omega_meas(step_cnt) = theta_dot + 0.5*(beta(step_cnt)+...
        beta(step_cnt-1))+randn*sqrt(sigma_nu^2/dt + (1/12)*sigma_u^2*dt);
    
    % propogate observer state prior to taking in measurments
    theta_hat_pri = theta_hat_post(step_cnt-1) +...
        omega_meas(step_cnt)*dt - beta_hat_post(step_cnt-1)*dt; 
    %     ^ intersting omega_meas is the current step...
    beta_hat_pri = beta_hat_post(step_cnt-1);
    
    % get e_meas from measurments
    e_meas = theta_meas(step_cnt) - theta_hat_pri;
    
    % update observer states after measurments
    theta_hat_post(step_cnt) = theta_hat_pri + K_theta*e_meas;
    beta_hat_post(step_cnt) = beta_hat_pri + K_beta*e_meas;
    
    e_truth(step_cnt) = theta_hat_post(step_cnt) - theta(step_cnt);
end
    
figure;
plot(t,e_truth);
title('Angular Error (Estimate - Truth) [rad]');
figure;
plot(t,theta_hat_post);

