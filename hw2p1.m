
% LM KF Class
% HW2 Problem 1
% 2018-11-10
% Implement Kalman Filter for 1-D tracking problem

%% Part a
clear all;
close all;

% sim parameters
dt = 1.0;
t_start = 0;
t_end = 300;
time = t_start:dt:t_end;

% truth parameters
Phi = [1 dt; 0 1];
H = [1 0];
sigma_nu = 0.1;

% kalman filter parameter
Q = [0 0; 0 0.001^2]; 
Gamma = eye(2);
R = 0.1^2;
Phi_kf = Phi; % setting state transition matrix Phi_kf same as actual truth, but might differ in reality
H_kf = H; % ditto ^


%% Propogate
for k = 1:length(time)
    if k == 1 % initialize state
        x_hat_prior(:,k) = [0.1 0.9]'; % x_hat_initial
        P_prior(:,:,k) = [0.5^2, 0; 0, 0.5^2]; % P_initial
        x(:,k) = [0 1]'; % truth initial state
    end
    
    % calculate K gain
    K(:,k) = P_prior(:,:,k)*H_kf'*inv(H_kf*P_prior(:,:,k)*H_kf' + R);
    
    % take in measurment
    y(:,k) = H*x(:,k) + sigma_nu*randn;
    % update state and covariance
    e_meas(:,k) = y(:,k) - H_kf*x_hat_prior(:,k);
    x_hat_post(:,k) = x_hat_prior(:,k) + K(:,k)*e_meas(:,k);
    P_post(:,:,k) = (eye(2) - K(:,k)*H_kf)*P_prior(:,:,k);
    
    e_states(:,k) = x(:,k) - x_hat_post(:,k); % record for plotting
    
    % propogate estimates forward in time
    x_hat_prior(:,k+1) = Phi_kf*x_hat_post(:,k);
    P_prior(:,:,k+1) = Phi_kf*P_post(:,:,k)*Phi_kf' + Gamma*Q*Gamma';
    one_sigma_pos(k) = sqrt(P_post(1,1,k));
    one_sigma_vel(k) = sqrt(P_post(2,2,k));
    
    % propogate truth
    x(:,k+1) = Phi*x(:,k);
end
    
%% Plot
figure()
ax1 = subplot(2,1,1);
plot(time,K(1,:));
ylabel('Position Gain')
ax2 = subplot(2,1,2);
plot(time,K(2,:));
ylabel('Velocity Gain')
linkaxes([ax1,ax2],'x')

figure()
ax1 = subplot(2,1,1);
plot(time,e_states(1,:));
hold on;
plot(time,3*one_sigma_pos,'r');
plot(time,-3*one_sigma_pos,'r');
ylim([-0.2 0.2]);
ylabel('Position Error');
legend({'position','3-\sigma'});
ax2 = subplot(2,1,2);
plot(time,e_states(2,:));
hold on;
plot(time,3*one_sigma_vel,'r');
plot(time,-3*one_sigma_vel,'r');
legend('velocity','3-\sigma');
ylabel('Velocity Error')
ylim([-0.02 0.02]);
linkaxes([ax1,ax2],'x')

figure()
ax1 = subplot(2,1,1);
plot(time,x(1,1:(end-1)));
hold on;
plot(time,x_hat_post(1,:),'r');
ylabel('Position(m)');
legend({'Position True (m)','Position Estimate (m)'});
ax2 = subplot(2,1,2);
plot(time,x(2,1:(end-1)));
hold on;
plot(time,x_hat_post(2,:),'r');
legend({'Velocity True (m/s)','Velocity Estimate (m/s)'});
ylabel('Velocity (m/s)')
linkaxes([ax1,ax2],'x')

%% Part b, truth state matrix modified
clear all;


% sim parameters
dt = 0.01;
t_start = 0;
t_end = 300;
time = t_start:dt:t_end;

% truth parameters
Phi = [1 dt; 0 0.95];
H = [1 0];
sigma_nu = 0.1;

% kalman filter parameter
Q = [0 0; 0 0.001^2]; 
Gamma = eye(2);
R = 0.1^2;
Phi_kf = Phi; % setting state transition matrix Phi_kf same as actual truth, but might differ in reality
H_kf = H; % ditto ^


%% Propogate
for k = 1:length(time)
    if k == 1 % initialize state
        x_hat_prior(:,k) = [0.1 0.9]'; % x_hat_initial
        P_prior(:,:,k) = 0.5^2*eye(2); % P_initial
        x(:,k) = [0 1]'; % truth initial state
    end
    
    % calculate K gain
    K(:,k) = P_prior(:,:,k)*H_kf'*inv(H_kf*P_prior(:,:,k)*H_kf' + R);
    
    % take in measurment
    y(:,k) = H*x(:,k) + sigma_nu*randn;
    % update state and covariance
    e_meas(:,k) = y(:,k) - H_kf*x_hat_prior(:,k);
    x_hat_post(:,k) = x_hat_prior(:,k) + K(:,k)*e_meas(:,k);
    P_post(:,:,k) = (eye(2) - K(:,k)*H_kf)*P_prior(:,:,k);
    
    e_states(:,k) = x(:,k) - x_hat_post(:,k); % record for plotting
    
    % propogate estimates forward in time
    x_hat_prior(:,k+1) = Phi_kf*x_hat_post(:,k);
    P_prior(:,:,k+1) = Phi_kf*P_post(:,:,k)*Phi_kf' + Gamma*Q*Gamma';
    one_sigma_pos(k) = sqrt(P_post(1,1,k));
    one_sigma_vel(k) = sqrt(P_post(2,2,k));
    
    % propogate truth
    x(:,k+1) = Phi*x(:,k);
end
    
%% Plot
figure()
ax1 = subplot(2,1,1);
plot(time,K(1,:));
ylabel('Position Gain')
% ylim([-0.5 0.5]);
ax2 = subplot(2,1,2);
plot(time,K(2,:));
ylabel('Velocity Gain')
% ylim([-0.01 0.05]);
linkaxes([ax1,ax2],'x')

figure()
ax1 = subplot(2,1,1);
plot(time,e_states(1,:));
hold on;
plot(time,3*one_sigma_pos,'r');
plot(time,-3*one_sigma_pos,'r');
ylim([-0.2 0.2]);
ylabel('Position Error');
legend({'position','3-\sigma'});
ax2 = subplot(2,1,2);
plot(time,e_states(2,:));
hold on;
plot(time,3*one_sigma_vel,'r');
plot(time,-3*one_sigma_vel,'r');
legend('velocity','3-\sigma');
ylabel('Velocity Error')
linkaxes([ax1,ax2],'x')

figure()
ax1 = subplot(2,1,1);
plot(time,x(1,1:(end-1)));
hold on;
plot(time,x_hat_post(1,:),'r');
ylabel('Position(m)');
legend({'Position True (m)','Position Estimate (m)'});
ax2 = subplot(2,1,2);
plot(time,x(2,1:(end-1)));
hold on;
plot(time,x_hat_post(2,:),'r');
legend({'Velocity True (m/s)','Velocity Estimate (m/s)'});
ylabel('Velocity (m/s)')
linkaxes([ax1,ax2],'x')


