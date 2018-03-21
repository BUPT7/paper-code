close all;
clc;
clear all;
clear global;
% % Install the SeDuMi
% path(path,'.\cvx-w64\cvx\sedumi');
% install_sedumi;
% % Install the Yalmip
% path(path,'.\');
% addpath(genpath('.\yalmip'));
%% Install the CVX
path(path,'.\cvx');
cvx_setup;
clc;
path(path,'.\');
fprintf('********************************************************************\n');
fprintf('cxy".\n');
fprintf('Communications and Network Center, School of Electronic Engineering.\n');
tic;
%% Init channel
N = 1;                                        % Numbers of Source Node's antennas
M = 1;                                        % Numbers of Destination Node's antennas
TTIs = 10;                                    % The number of channel pairs of this simulation
str = strcat('./channels/H_',num2str(N),'_',num2str(M),'_',num2str(TTIs),'.mat');
load(str);
clear -regexp [^H ^N ^M ^TTIs];

%% init parameters
dLink = 0.1:0.2:2.5;
iter_num =30;
beta1 = 0.5;
beta2 = 0.5;
schemes = [ 1,2 ];

%% power constraints
P1_max = 30;
P2_max = 30;
e1 = 0.5;
e2 = 0.5;

%% init accumulator
number = 2;
Average_rate = zeros(length(dLink),number,iter_num+1);
Average_p1 = zeros(length(dLink),number,iter_num+1);
Average_p2 = zeros(length(dLink),number,iter_num+1);

%% Simulation starts
% matlabpool('open');    
for d_L = 1:length( dLink )
dL = dLink( d_L );
fprintf('The %d th calculation starts at dL = %3.2f ... \n',d_L,dL);
[ rate,p1,p2] = monte_carlo( TTIs,H,dL,N,M,iter_num,P1_max,P2_max,e1,e2,beta1,beta2,number,schemes );
Average_rate( d_L,:,:) = rate / TTIs;
Average_p1( d_L,:,:) = p1 / TTIs;
Average_p2( d_L,:,:) = p2 / TTIs;
end
% matlabpool('close');

results_folder = './results';
results_file = fullfile(results_folder,...
    sprintf('Chen_GAME_dL=14_TTI=10_p=30_sigma=1_p1=1_p2=0.5.mat'));
save(results_file);
fprintf('Saving results to ./results/...\n');
toc;