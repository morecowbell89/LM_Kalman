% LM KF Class
% HW2 Problem 2
% 2018-11-10
% Implement Kalman Filter for Gyro

%% 
clear all;
close all;

% sim parameters
dt = 0.01;
t_start = 0;
t_end = 1000;
time = t_start:dt:t_end;

% truth parameters
sigma_nu = 3e-7; % Angle random walk (rad/sec^0.5)
sigma_u = 3e-10; % rate random walk (rad/sec^(3/2))
sigma_S = 17e-6; % noise on measured angle
omega_in = 1e-7; % constant rate truth (rad/sec)
Phi = [1 0; 0 1]; % calculate state transition model; states are [theta beta]
Lambda = [dt 0]'; % used for inputing omega*dt into theta
Gamma = [0 sigma_u*sqrt(dt)]';
H = [1 0];

% kalman filter parameter
Q = [(sigma_nu^2*dt + (1/3)*sigma_u^2*dt^3) -0.5*sigma_u^2*dt; -0.5*sigma_u^2*dt sigma_u^2*dt]; 
Gamma_kf = [dt 0; 0 0];
R = sigma_S^2;
Phi_kf = [1 -dt; 0 1]; 
H_kf = H; % ditto ^

%% Propogate
for k = 1:length(time)
    if k == 1 % initialize state
        x_hat_prior(:,k) = [0 0]'; % x_hat_initial
        P_prior(:,:,k) = (1e-4)*eye(2); % P_initial
        x(:,k) = [0 1e-6]'; % truth initial state [theta beta]'
    end
    omega(k) = omega_in; % flexible for time varying. constant though
    omega(k+1) = omega_in; % flexible for time varying. constant
    
    % propogate truth
    x(:,k+1) = Phi*x(:,k) + Lambda*omega(k+1) + Gamma*randn;
    % calculate K gain
    K(:,k) = P_prior(:,:,k)*H_kf'*inv(H_kf*P_prior(:,:,k)*H_kf' + R);
    
    % take in measurment
    y(:,k) = H*x(:,k) + sigma_S*randn; % theta
    omega_meas = omega(k+1) + 0.5*(x(2,k+1)+x(2,k)) + randn*sqrt(sigma_nu^2/dt + (sigma_u^2*dt)/12);
    % update state and covariance
    e_meas(:,k) = y(:,k) - H_kf*x_hat_prior(:,k);
    x_hat_post(:,k) = x_hat_prior(:,k) + K(:,k)*e_meas(:,k);
    P_post(:,:,k) = (eye(2) - K(:,k)*H_kf)*P_prior(:,:,k);
    
    e_states(:,k) = x(:,k) - x_hat_post(:,k); % record for plotting
    
    % propogate estimates forward in time
    x_hat_prior(:,k+1) = Phi_kf*x_hat_post(:,k) + Lambda*omega_meas;
    P_prior(:,:,k+1) = Phi_kf*P_post(:,:,k)*Phi_kf' + Gamma_kf*Q*Gamma_kf';
    one_sigma_pos(k) = sqrt(P_post(1,1,k));
    one_sigma_vel(k) = sqrt(P_post(2,2,k));
    
    
end

%% Plot
figure()
ax1 = subplot(2,1,1);
plot(time,K(1,:));
ylabel('Theta Gain')
ax2 = subplot(2,1,2);
plot(time,K(2,:));
ylabel('Beta Gain')
linkaxes([ax1,ax2],'x')

figure()
ax1 = subplot(2,1,1);
plot(time,e_states(1,:));
hold on;
plot(time,3*one_sigma_pos,'r');
plot(time,-3*one_sigma_pos,'r');
% ylim([-0.2 0.2]);
ylabel('\theta Error');
legend({'\theta','3-\sigma'});
ax2 = subplot(2,1,2);
plot(time,e_states(2,:));
hold on;
plot(time,3*one_sigma_vel,'r');
plot(time,-3*one_sigma_vel,'r');
legend('\beta','3-\sigma');
ylabel('\beta Error')
% ylim([-0.02 0.02]);
linkaxes([ax1,ax2],'x')

figure()
ax1 = subplot(2,1,1);
plot(time,x(1,1:(end-1)));
hold on;
plot(time,x_hat_post(1,:),'r');
ylabel('\theta (rad)');
legend({'\theta True (m)','\theta Estimate (m)'});
ax2 = subplot(2,1,2);
plot(time,x(2,1:(end-1)));
hold on;
plot(time,x_hat_post(2,:),'r');
legend({'\beta True (m/s)','\beta Estimate (m/s)'});
ylabel('\beta (rad/s)')
linkaxes([ax1,ax2],'x')
