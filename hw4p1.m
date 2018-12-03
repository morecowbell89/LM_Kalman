% Kalman Filter Class
% HW 4, Question 1
% 2018-12-04

% 3-axis dynamic model
clear all;close all;

% Dynamics parameters
Inertia = [100, 5, 8;
            5, 120, 9;
            8,9,150]; % Kg-m^2 

invInertia = inv(Inertia);
        
w0 = [5e-3, 6e-3, 7e-3]'; % rad/sec

qi2b0 = [0,1,0,0]';

Torque = zeros(3,1);  

% Gyro Parameters
sigmaV = 3e-7;      % Angle random walk (rad/sec^0.5)
sigmaU = 3e-10;     % Rate random walk (rad/sec^(3/2) )
beta = (1e-6)*[0.5,1,-1.5]';    % initial bias x,y,z (rad/sec)

% Inertial Measurement noise
sigmaECI = 0.02;

% Initial truth state - quaternion and rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [qi2b0;w0];
x = x0;

v1(:,1) = q2dcm(x(1:4,1))*[1,0,0]'+[randn;randn;randn]*sigmaECI;
v2(:,1) = q2dcm(x(1:4,1))*[0,1,0]'+[randn;randn;randn]*sigmaECI;

dt = 0.1;
time = 0:dt:300;

% Initial estimate state - quaternion and rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_non_norm = [0.02 1 0.02 0]';
q_hat = q_non_norm./norm(q_non_norm);
b_hat = [0 0 0]';
delta_x_hat = zeros(6,1); % [a1 a2 a3 b1 b2 b3]'
P_kf = 1e-4*eye(6); % Covariance
% Process Noise
Q_kf= [(sigmaV^2*dt + (1/3)*sigmaU*dt^3)*eye(3) ((1/2)*sigmaU^2*dt^2)*eye(3);
                 ((1/2)*sigmaU^2*dt^2)*eye(3)         (sigmaU^2*dt)*eye(3)];
H_kf = [eye(3) zeros(3); eye(3) zeros(3)]; 
R_kf = sigmaECI^2*eye(6);
Gamma_kf = [-eye(3) zeros(3); zeros(3) eye(3)];

% capture intial estimates to spool
q_hat_spool(:,1) = q_hat;
b_hat_spool(:,1) = b_hat;
one_sigma_x(1) = sqrt(P_kf(1,1));
z_spool_cnt = 1;
% Start Sim/Kalman Filter
for k = 1:length(time)
    % Truth Sim Section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Simulate rigid body truth model
    % Torque = I*wdot + omega x Iw
    % wdot = inv(I)*(torque - cross(omega, Iw))
    % Calculate body momentum
    H(:,k) = Inertia*x(5:7,k);
    H_ECI(:,k) = q2dcm(qinv(x(1:4,k)))*H(:,k);
    % Compute angular acceleration
    wdot = inv(Inertia)*(Torque - cross(x(5:7,k),H(:,k)));
    % Compute quaternion derivative
    qdot = qmult(0.5*[x(5:7,k);0],x(1:4,k));
    
    % Formulate state derivative and propagate
    xDot = [qdot;wdot];
    x(:,k+1) = x(:,k) + xDot*dt; % Propagate state
    
    % Generate gyro measurements
    % Propagate true gyro bias state
    beta(:,k+1) = beta(:,k) + [randn;randn;randn]*sigmaU*sqrt(dt);
    
    % Calculate angular rate measurement from gyro
    omegaTilde(:,k) = x(5:7,k+1) + 0.5*(beta(:,k+1)+beta(:,k)) + randn*sqrt( sigmaV^2/dt + (1/12)*sigmaU^2*dt );    % changed this from bill's
  
    % At 1 Hz, compute measurement vectors in body frame
    if mod(time(k),1)==0
        v1(:,k+1) = q2dcm(x(1:4,k+1))*[1,0,0]'+[randn;randn;randn]*sigmaECI;
        v2(:,k+1) = q2dcm(x(1:4,k+1))*[0,1,0]'+[randn;randn;randn]*sigmaECI;
    else
        v1(:,k+1) = v1(:,k);
        v2(:,k+1) = v2(:,k);
    end
    
    % Kalman Section
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1) Kalman Gain
    %%%%%%%%%%%%%%%%%%%%%%%%%
    K_kf = P_kf*H_kf'*inv(H_kf*P_kf*H_kf' + R_kf);
    
    % 2) Update with measurments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(time(k),1)==0
        z(1:3,1) = v1(:,k) - q2dcm(q_hat)*[1,0,0]';
        z(4:6,1) = v2(:,k) - q2dcm(q_hat)*[0,1,0]';
        
        delta_x_hat = K_kf*z;
        q_hat = qmult([0.5*delta_x_hat(1:3);1],q_hat);
        b_hat = delta_x_hat(4:6) + b_hat;
        P_kf = (eye(6)-K_kf*H_kf)*P_kf;
        
        z_spool(:,z_spool_cnt) = z(:,1);
        z_spool_cnt = z_spool_cnt + 1;
    end
    
    % Capture to spool
    if k > 1
        q_hat_spool(:,k) = q_hat;
        b_hat_spool(:,k) = b_hat;
        one_sigma_x(k) = sqrt(P_kf(1,1));
        
    end
    
    % 3) Propagate to next time step
    %%%%%%%%%%%%%%%%%%%%%%%
    % compensate the gyro data with estimated bias
    delta_theta_hat = omegaTilde(:,k)*dt - b_hat*dt;
    q_hat = qmult([0.5*delta_theta_hat;(1-(1/8)*dot(delta_theta_hat,delta_theta_hat))],q_hat);
    temp =  [0 -delta_theta_hat(3) delta_theta_hat(2); 
             delta_theta_hat(3) 0 -delta_theta_hat(1);
             -delta_theta_hat(2) delta_theta_hat(1) 0];
    Phi_one = eye(3) + temp;
    Phi_kf = [Phi_one -eye(3)*dt; zeros(3) eye(3)];
    P_kf = Phi_kf*P_kf*Phi_kf' + Gamma_kf*Q_kf*Gamma_kf';    
    
    
    
end
subplot(2,1,1)
plot(time,q_hat_spool)
subplot(2,1,2)
plot(time,x(1:4,1:end-1))

figure()
plot(time,one_sigma_x);

figure
plot(time,omegaTilde(:,1:end))
title('Gyro Measurement');
xlabel('Time (s)');
ylabel('Rad/sec');
legend('x','y','z');