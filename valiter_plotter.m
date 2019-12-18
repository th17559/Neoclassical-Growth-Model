%% This script plots the value function and the policy function 

clc;
clear all;
%% Read file
% w. interpolation
fileID = fopen('valiter_cspline');
X=fscanf(fileID,'%f %f %f %f %f',[5 Inf]);
X=X';

k_cs = X(:,1);
v1_cs = X(:,2);
v2_cs = X(:,3);
g1_cs = X(:,4);
g2_cs = X(:,5);

%% Plot value and policy functions

figure(1)
% plot value functions
plot(k_cs,v1_cs,k_cs,v2_cs)
legend('z=1.05','z=0.84')

figure(2)
% plot policy functions
plot(k_cs,g1_cs,k_cs,g2_cs)
hold on
plot(k_cs,k_cs,'--')
legend('z=1.05','z=0.84')

%% Plot derivative of value and policy functions

% define anonymous functions
v1 = @(k) interp1(k_cs,v1_cs,k,'spline') ;
v2 = @(k) interp1(k_cs,v2_cs,k,'spline') ;
g1 = @(k) interp1(k_cs,g1_cs,k,'spline') ;
g2 = @(k) interp1(k_cs,g2_cs,k,'spline') ;

% approximate derivative
h = 0.001;
klb = 1;    % klb set to 1 since derivative of value function for klb=1e-5 is a very large number
kub = 45;
K = klb:h:kub;

v1p = diff(v1(K))/h ;
v2p = diff(v2(K))/h ;

figure(3)
% plot derivative of value functions
plot(K(1:length(v1p)),v1p, K(1:length(v2p)),v2p)
legend('z=1.05','z=0.84')

g1p = diff(g1(K))/h ;
g2p = diff(g2(K))/h ;

figure(4)
% plot derivative of value functions
plot(K(1:length(g1p)),g1p, K(1:length(g2p)),g2p)
legend('z=1.05','z=0.84')