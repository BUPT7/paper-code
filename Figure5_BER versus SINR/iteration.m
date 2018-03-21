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
aLink =1:1:10;%接收功率阈值
%aL=0.1;
iter_num =10;
beta1 = 0.5; 
beta2 = 0.5;
beta3 = 0.5;
beta4 = 0.5;
schemes = [1,0,0,0];

%% power constraints
P1_max = 1;
P2_max = 1;
P3_max = 1;
P4_max = 1;
dL = 1;
e = 0.5;

%% init accumulator
number = 4;
Average_rate = zeros(length(aLink),number,iter_num+1);
Average_p1 = zeros(length(aLink),number,iter_num+1);
Average_p2 = zeros(length(aLink),number,iter_num+1);
Average_p3 = zeros(length(aLink),number,iter_num+1);
Average_p4 = zeros(length(aLink),number,iter_num+1);
%Average_p = zeros(length(aLink),number,iter_num+1);

%% Simulation starts
% matlabpool('open');    
for a_L = 1:length( aLink )
aL = aLink( a_L );
fprintf('The %d th calculation starts at eL = %3.2f ... \n',a_L,aL);
[ rate,p1,p2,p3,p4,p] = monte_carlo( TTIs,H,dL,N,M,iter_num,P1_max,P2_max,P3_max,P4_max,aL,beta1,beta2,beta3,beta4,number,schemes,e );
Average_rate( a_L,:,:) = rate / TTIs;
%Average_p( a_L,:,:) = p / TTIs;
Average_p1( a_L,:,:) = p1 / TTIs;
Average_p2( a_L,:,:) = p2 / TTIs;
Average_p3( a_L,:,:) = p3 / TTIs;
Average_p4( a_L,:,:) = p4 / TTIs;
end
% matlabpool('close');

results_folder = './results';
results_file = fullfile(results_folder,...
    sprintf('Chen_GAME_dL=14_TTI=10_p=30_sigma=1_p1=1_p2=0.5.mat'));
save(results_file);
fprintf('Saving results to ./results/...\n');
toc;