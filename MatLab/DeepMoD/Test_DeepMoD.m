%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%         This Matlab File uses the S^3 method described in           %%%
%%%     Ye Yuan, "Discovery of Partial Differential Equations from      %%%
%%%                        Spatiotemporal Data"                         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           S^3d is applied to analyze the following                  %%%
%%%                     heat transfert equation:                        %%%
%%%          u{t} = Dc*u{xx} - UL*(u{xx}-u_amb) + Q(x,t)                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
addpath Functions
%% Generate data

% load Donnee_sim_Null_Flow.mat
% load Time & Space.mat
data = load('Donnee_sim_Null_Flow.mat');
u = data.data; 
dt = (data.Time(2) - data.Time(1))*3600;
dx = data.Layers(2) - data.Layers(1);
t = data.Time*3600;
x = data.Layers;
%% Estimate derivatives

% %Method: polynomial interpolation
% t_th = 1; %order of time derivative
% x_th = 3; %maximum order of spatial derivative
%        
% parameter.deg = 3;           %degree of polynomial to use 
% parameter.x_num_to_fit = 3; %number of points to use in polynomial interpolation for spatial derivatives
% parameter.t_num_to_fit = 100; %number of points to use in polynomial interpolation for time derivative
% 
% [derivative] = make_input_poly(u,x,x_th,parameter);
% y0 =  make_y_poly(u,t,t_th,parameter );
% 
% input = derivative.derivative;
% input_name = derivative.name;
% 
%  %Add a reshaped version of Q_mat to input matrix, before making our theta
% [r,v] = size(data.Q_mat);
% num = (r-2*parameter.x_num_to_fit)*(v-2*parameter.t_num_to_fit);
% Q_vec  = reshape(data.Q_mat(1+parameter.x_num_to_fit:end-parameter.x_num_to_fit, 1+parameter.t_num_to_fit:end-parameter.t_num_to_fit), num, 1);
% 
% input = [input Q_vec];
% input_name = [input_name 'Q'];
% clearvars r v num

%% Estimate derivatives

%Method: Finite difference
t_th = 1; %order of time derivative
x_th = 2; %maximum order of spatial derivative

[derivative] = make_input_fd( u,dx,x_th);        
y0 = make_y_fd( u,dt,t_th );

input = derivative.derivative;
input_name = derivative.name;

 %Add a reshaped version of Q_mat to input matrix, before making our theta

Q_vec  = reshape(data.Q_mat, [], 1);
input = [input Q_vec];
input_name = [input_name 'Q'];

%% Builds a dictionary matrix

polyorder = 1; %maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);

 %Add a reshaped version of Q_mat to input matrix, before making our theta

% Q_vec  = reshape(data.Q_mat, [], 1);
% theta0 = [theta0 Q_vec];
% theta_name = [theta_name 'Q'];

%Restreignont le dictionnaire aux valeurs qu'on veux(et qu'on connais, c'est de la triche mais bon ........)
I = strcmp(theta_name, '1') | strcmp(theta_name, 'u') | strcmp(theta_name, 'u{xx}') | strcmp(theta_name, 'Q');
theta0 = theta0(:,I);
theta_name = theta_name(I);



%% Number of snapshots: 5000

index = randperm(6000,5000);
theta = theta0(index,:);
y = y0(index,:);

%% Compute Sparse regression: Sparse Bayesian Approach

%normalization
T = zeros(size(theta,1),size(theta,2));
Mreg = zeros(size(theta,2),1);
for i = 1:size(theta,2)
    Mreg(i,1) = 1/max(abs(theta(:,i)));
    T(:,i) = Mreg(i)*theta(:,i);    
end

lambda = 1e-2; %the regularization parameter 
MAXITER = 1000;   %maximum number of inner iterations
w = tac_reconstruction(y, T, lambda,MAXITER);
w = w(:,end).*Mreg; %identified parameters corresponding to the basis function in Theta

%% print result

threshold = 1e-8; %filter those cofficients in w being less than the threshold
y_name = 'u{t}';
fprintf('\n%s = ', y_name);
for i = 1:size(w,1)
    if abs(w(i))<threshold
        w(i) = 0;
    else
        if w(i)<0
            fprintf('%f%s', w(i),theta_name{i});
        else
            fprintf('+');
            fprintf('%f%s', w(i),theta_name{i});
        end
    end    
end
%% print error
% u{t} = Dc*u{xx} - UL*(u{xx}-u_amb) + Q(x,t)
% RMS Error
err = abs([(data.Tank.Dc*data.eps -w(3))/(data.Tank.Dc*data.eps), (-data.Tank.UL*10 -w(2))/(data.Tank.UL*10), (1 - w(4))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);