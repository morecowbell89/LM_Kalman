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
Phi = eye(6) + A*dt + A.^2*dt; % not sure if it should be A.^2 or A*A

% Process Noise
B = [0 ,0 ,0;
     0 ,0 ,0;
     0 ,0 ,0;
     1 ,0 ,0;
     0 ,1 ,0;
     0 ,0 ,1;];
W = sigma_proc^2*eye(3);
S = B*W*B';
Q = 0.5*dt*(S+Phi*S*Phi');

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
K = zeros(3,length(time));
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
         P_prior(:,:,k) = 5*sigma_proc*eye(6);
     end
     
     % Propogate truth states
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x(:,k+1) = Phi*x(:,k);
     
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
     y(:,k) = H*x(:,k) + [sigma_alpha*randn sigma_beta*randn sigma_r*randn]';
     
     % Calculate K gain
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Calculate H_kf, using estimated state
     r_kf = sqrt(x_hat_prior(1,k)^2+x_hat_prior(2,k)^2+x_hat_prior(3,k)^2);
     H(1,1) = x_hat_prior(3,k)/(x_hat_prior(1,k)^2 + x_hat_prior(3,k)^2);
     H(1,2) = 0;
     H(1,3) = -x_hat_prior(1,k)/(x_hat_prior(1,k)^2 + x_hat_prior(3,k)^2);
     H(2,1) = -(x_hat_prior(1,k)*x_hat_prior(2,k))/(r_kf^2*sqrt(x_hat_prior(1,k)^2+x_hat_prior(3,k)^2));
     H(2,2) = sqrt(x_hat_prior(1,k)^2 + x_hat_prior(3,k)^2)/r_kf^2;
     H(2,3) = -(x_hat_prior(2,k)*x_hat_prior(3,k))/(r_kf^2*sqrt(x_hat_prior(1,k)^2+x_hat_prior(3,k)^2));
     H(3,1) = x_hat_prior(1,k)/r_kf;
     H(3,2) = x_hat_prior(2,k)/r_kf;
     H(3,3) = x_hat_prior(3,k)/r_kf;
     % Calculate K
     K(:,k) = P_prior(:,:,k)*H_kf'*inv(H_kf*P_prior(:,:,k)*H_kf' + R);
     
     % Calculate innovation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     e_meas(:,k) = y(:,k) - H_kf*x_hat_prior(:,k);
     
     % Apply corrections
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % update state posteriori
     x_hat_post(:,k) = x_hat_prior(:,k) + K(:,k)*e_meas(:,k);

