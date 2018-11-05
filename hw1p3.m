clear all;
close all;
clc;

req_deg = 1; % [deg] pointing requirment
req_rad = req_deg*(pi/180);

% gyro statistics
sigma_nu = 3e-7; %[rad/sec^0.5]
sigma_u = 3e-10; %[rad/sec^1.5]

% 1 sigma equation: sqrt((1/3)*sigma_u^2*t^3+sigma_nu^2*t)
% retain confidence within 3sigma of point requirement as:
% 3*sqrt((1/3)*sigma_u^2*t^3+sigma_nu^2*t) = req_rad
% solve for zeros with:
% 3*sqrt((1/3)*sigma_u^2*t^3+sigma_nu^2*t) - req_rad = 0

confidence_calc = @(t) 3*sqrt((1/3)*sigma_u^2*t.^3+sigma_nu^2*t) - req_rad;

% seed fzero
t0 = 1e5;

t_confidence = fzero(confidence_calc,t0); % [seconds]

t_confidence_hours = t_confidence/3600;

fprintf('Can propogate from gyro for %4.2f hours while being within 3-sigma of +/- 1.0deg pointing requirement \n',t_confidence_hours);