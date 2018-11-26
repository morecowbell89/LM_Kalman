% LM KF Class
% HW3 Problem 2
% 2018-11-25
% Implement Satellite Relative Navigation Tracking

%% Setting Parameters/Allocating vectors
clear all;
close all;

% sim time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1; % in order for state equations to work dt needs to be in units of [sec]
t_start = 0;
t_end = 6000;
time = t_start:dt:t_end;

% truth parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orbit_period = 1.5*3600; % 1.5 hours*3600 seconds [sec]
n = 2*pi/orbit_period; % [rad/sec]
sigma_alpha = 1*(pi/180); % measurment noise 1-sigma in alpha [rad]
sigma_beta = 1*(pi/180); % measurement noise 1-sigma in beta [rad]
sigma_r = 2*(1/100); % measurment noise in range 1-sigma [meter]
sigma_proc = 1e-4; % process noise for all 3-axis 1-sigma [m/s^2]

% State Transition Matrix
A = [0     ,0   ,0    ,1    ,0   ,0;
     0     ,0   ,0    ,0    ,1   ,0;
     0     ,0   ,0    ,0    ,0   ,1;
     3*n^2 ,0   ,0    ,0    ,2*n ,0;
     0     ,0   ,0    ,-2*n ,0   ,0;
     0     ,0   ,-n^2 ,0    ,0   ,0;];
Phi = eye(6) + A*dt + A^2*(dt^2/2); % not sure if it should be A.^2 or A*A

% Process Noise
B = [0 ,0 ,0;
     0 ,0 ,0;
     0 ,0 ,0;
     1 ,0 ,0;
     0 ,1 ,0;
     0 ,0 ,1;];
W = sigma_proc^2*eye(3);
S = B*W*B';

% Measurment noise
R = [sigma_alpha^2 ,0            ,0         ;
     0             ,sigma_beta^2 ,0         ;
     0             ,0            ,sigma_r^2 ;];
 
% preallocating vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(6,length(time)+1);
x_hat_prior = zeros(6,length(time)+1);
P_prior = zeros(6,6,length(time));
H = zeros(3,6);
H_kf = zeros(3,6);
y = zeros(3,length(time));
K = zeros(6,3,length(time));
e_meas = zeros(3, length(time));
 %% Run Simulation
 for k = 1:length(time)
     if k == 1
         % Initialize states
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % initialize truth state
         x(:,k) = [1 0.1 0.1 0 0 0]'; % truth state [m,m,m,m/s,m/s,m/s]
         % initialize state estimate
         x_hat_prior(:,k) = x(:,k); % set initial estimate equal to truth
         % initialize covariance estimate
         P_prior(:,:,k) = 20*sigma_proc*eye(6);
     end
     
     % Propogate truth states
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x(:,k+1) = Phi*x(:,k)+B*randn(3,1)*sigma_proc;
     
     % Take in measurments
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % compute truth angles measured based off of truth staetes
     r = sqrt(x(1,k)^2+x(2,k)^2+x(3,k)^2);
     H(1,1) = x(3,k)/(x(1,k)^2 + x(3,k)^2);
     H(1,2) = 0;
     H(1,3) = -x(1,k)/(x(1,k)^2 + x(3,k)^2);
     H(2,1) = -(x(1,k)*x(2,k))/(r^2*sqrt(x(1,k)^2+x(3,k)^2));
     H(2,2) = sqrt(x(1,k)^2 + x(3,k)^2)/r^2;
     H(2,3) = -(x(2,k)*x(3,k))/(r^2*sqrt(x(1,k)^2+x(3,k)^2));
     H(3,1) = x(1,k)/r;
     H(3,2) = x(2,k)/r;
     H(3,3) = x(3,k)/r;
     % other elements of H equal 0 always since the elements are partial 
     % derivatives w.r.t. velocity components, and observations independent
     % of velocity
     
%      [atan(-x(3,k)./x(1,k)), asin(x(2,k)/norm(x(:,k))), norm(x(:,k)) ]'
% + [sigma_alpha*randn sigma_beta*randn sigma_r*randn]'
     y(:,k) =  H*x(:,k);
     
     % Calculate K gain
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Calculate H_kf, using estimated state
     r_kf =  norm(x_hat_prior(:,k));
     H_kf(1,1) = x_hat_prior(3,k)/(x_hat_prior(1,k)^2 + x_hat_prior(3,k)^2);
     H_kf(1,2) = 0;
     H_kf(1,3) = -x_hat_prior(1,k)/(x_hat_prior(1,k)^2 + x_hat_prior(3,k)^2);
     H_kf(2,1) = -(x_hat_prior(1,k)*x_hat_prior(2,k))/(r_kf^2*sqrt(x_hat_prior(1,k)^2+x_hat_prior(3,k)^2));
     H_kf(2,2) = sqrt(x_hat_prior(1,k)^2 + x_hat_prior(3,k)^2)/r_kf^2;
     H_kf(2,3) = -(x_hat_prior(2,k)*x_hat_prior(3,k))/(r_kf^2*sqrt(x_hat_prior(1,k)^2+x_hat_prior(3,k)^2));
     H_kf(3,1) = x_hat_prior(1,k)/r_kf;
     H_kf(3,2) = x_hat_prior(2,k)/r_kf;
     H_kf(3,3) = x_hat_prior(3,k)/r_kf;
     % Calculate K
     K(:,:,k) = P_prior(:,:,k)*H_kf'*inv(H_kf*P_prior(:,:,k)*H_kf' + R);
     
     % Calculate innovation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     e_meas(:,k) = y(:,k) - H_kf*x_hat_prior(:,k);
     
     % Apply corrections
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % update state posteriori
     x_hat_post(:,k) = x_hat_prior(:,k) + K(:,:,k)*e_meas(:,k);
     % update covariance posteriori
     P_post(:,:,k) = (eye(6) - K(:,:,k)*H_kf)*P_prior(:,:,k);
     % record error for plotting
     e_states(:,k) = x(:,k) - x_hat_post(:,k);
     
     % Propogate state and covariance to next time step (k+1)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x_hat_prior(:,k+1) = Phi*x_hat_post(:,k);
     P_prior(:,:,k+1) = Phi*(P_post(:,:,k)+0.5*dt*S)*Phi' + 0.5*dt*S;
     
     % calculate one sigma values for x,y,z,x',y',z'
     one_sigma_x(k) = sqrt(P_post(1,1,k));
     one_sigma_y(k) = sqrt(P_post(2,2,k));
     one_sigma_z(k) = sqrt(P_post(3,3,k));
     one_sigma_x_prime(k) = sqrt(P_post(4,4,k));
     one_sigma_y_prime(k) = sqrt(P_post(5,5,k));
     one_sigma_z_prime(k) = sqrt(P_post(6,6,k));
 end
 
 %% Plots
 %%%%%%%%%%%%%%%%%%%%%%
figure()
ax = subplot(2,1,1);
plot(time,x(1:3,1:(end-1)));
hold on;
ax.ColorOrderIndex = 1;
% plot(time,x_hat_post(1:3,:),'--');
ylabel('Position (m)');
xlabel('Time (sec)');
legend('x','y','z','xHat','yHat','zHat');
ax = subplot(2,1,2);
plot(time,x(4:6,1:(end-1)));
hold on;
ax.ColorOrderIndex = 1;
% plot(time,x_hat_post(4:6,:),'--');
ylabel('Velocity (m/s)');
xlabel('Time (sec)');
legend('x','y','z','xHat','yHat','zHat');

figure()
subplot(2,1,1)
plot(time,y(1:2,:)*180/pi)
ylabel('Angles (deg)');
legend('alpha','beta')
subplot(2,1,2)
plot(time,y(3,:));
ylabel('Range (m)');

figure()
subplot(2,1,1)
plot(time,e_meas(1:2,:)*180/pi)
ylabel('Angles (deg)');
legend('alpha','beta')
subplot(2,1,2)
plot(time,e_meas(3,:));
ylabel('Range (m)');

 
 

