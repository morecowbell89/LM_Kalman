%% Low Gain K
clear all;
close all;
clc;

dt = 1; % [sec]
Phi_hat = [1 dt; 0 1];
Phi = [1 dt; 0 1];
H = [1 0];
sigma_nu = 0.1;
K = [0.05; 0.03]; % low

cnt = 1;
x_hat_pri = [0.1 0.9]'; 
x(:,cnt) = [0 1]';
y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);

t_end = 300;
for cnt = 2:dt:(t_end+1)
    x(:,cnt) = Phi*x(:,cnt-1);
    y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
    
    x_hat_pri = Phi_hat*x_hat_post(:,cnt-1);
    x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
    est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);
end

t = (0:dt:t_end);

subplot(2,1,1)
plot(t,x(1,:))
hold on;
plot(t,x_hat_post(1,:));
legend('x','xhat');
ylabel('Position(m)');
xlabel('Time(sec)');
subplot(2,1,2)
plot(t,x(2,:))
hold on;
plot(t,x_hat_post(2,:));
legend('x','xhat');
ylabel('Velocity(m/sec)');
xlabel('Time(sec)');

%%  Medium Gain K
clear all;
close all;
clc;

dt = 1; % [sec]
Phi_hat = [1 dt; 0 1];
Phi = [1 dt; 0 1];
H = [1 0];
sigma_nu = 0.1;
K = [0.5; 0.3]; % low

cnt = 1;
x_hat_pri = [0.1 0.9]';
x(:,cnt) = [0 1]';
y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);

t_end = 300;
for cnt = 2:dt:(t_end+1)
    x(:,cnt) = Phi*x(:,cnt-1);
    y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
    
    x_hat_pri = Phi_hat*x_hat_post(:,cnt-1);
    x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
    est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);
end

t = (0:dt:t_end);

subplot(2,1,1)
plot(t,x(1,:))
hold on;
plot(t,x_hat_post(1,:));
legend('x','xhat');
ylabel('Position(m)');
xlabel('Time(sec)');
subplot(2,1,2)
plot(t,x(2,:))
hold on;
plot(t,x_hat_post(2,:));
legend('x','xhat');
ylabel('Velocity(m/sec)');
xlabel('Time(sec)');
%% High Gain K
clear all;
close all;
clc;

dt = 1; % [sec]
Phi_hat = [1 dt; 0 1];
Phi = [1 dt; 0 1];
H = [1 0];
sigma_nu = 0.1;
K = [1; 0.5]; % low

cnt = 1;
x_hat_pri = [0.1 0.9]';
x(:,cnt) = [0 1]';
y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);

t_end = 300;
for cnt = 2:dt:(t_end+1)
    x(:,cnt) = Phi*x(:,cnt-1);
    y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
    
    x_hat_pri = Phi_hat*x_hat_post(:,cnt-1);
    x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
    est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);
end

t = (0:dt:t_end);

subplot(2,1,1)
plot(t,x(1,:))
hold on;
plot(t,x_hat_post(1,:));
legend('x','xhat');
ylabel('Position(m)');
xlabel('Time(sec)');
subplot(2,1,2)
plot(t,x(2,:))
hold on;
plot(t,x_hat_post(2,:));
legend('x','xhat');
ylabel('Velocity(m/sec)');
xlabel('Time(sec)');
%% Bonus
clear all;
close all;
clc;

dt = 1; % [sec]
Phi_hat = [1 dt; 0 1]; % Phi_hat slightly different than Phi
Phi = [1 dt; 0 0.95];
H = [1 0];
sigma_nu = 0.1;
K = [1; 0.9]; % low

cnt = 1;
x_hat_pri = [0.5 0.95]';
x(:,cnt) = [0 1]';
y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);

t_end = 300;
for cnt = 2:dt:(t_end+1)
    x(:,cnt) = Phi*x(:,cnt-1);
    y(:,cnt) = H*x(:,cnt) + sigma_nu*randn;
    
    x_hat_pri = Phi_hat*x_hat_post(:,cnt-1);
    x_hat_post(:,cnt) = x_hat_pri + K*(y(:,cnt)-H*x_hat_pri);
    est_err(:,cnt) = x(:,cnt) - x_hat_post(:,cnt);
end

t = (0:dt:t_end);

subplot(2,1,1)
plot(t,x(1,:))
hold on;
plot(t,x_hat_post(1,:));
legend('x','xhat');
ylabel('Position(m)');
xlabel('Time(sec)');
subplot(2,1,2)
plot(t,x(2,:))
hold on;
plot(t,x_hat_post(2,:));
legend('x','xhat');
ylabel('Velocity(m/sec)');
xlabel('Time(sec)');
