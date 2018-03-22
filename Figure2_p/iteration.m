close all;
clc;
clear all;
clear global;
tic;

%% Init channel
N = 1;                                        % Numbers of Source Node's antennas
M = 1;                                        % Numbers of Destination Node's antennas
TTIs = 10;                                    % The number of channel pairs of this simulation
str = strcat('./channels/H_',num2str(N),'_',num2str(M),'_',num2str(TTIs),'.mat');
load(str);
clear -regexp [^H ^N ^M ^TTIs];

%% init parameters
dLink = 0.1:0.4:2.5;
iter_num = 20;
beta1 = 0.5;
beta2 = 0.5;
beta3 = 0.5;
beta4 = 0.5;
schemes = [ 1,2 ];

%% power constraints
P1_max = 30;
P2_max = 30;
P3_max = 30;
P4_max = 30;
e1 = 0.5;
e2 = 0.5;
e3 = 0.5;
e4 = 0.5;

%% init accumulator
number = 2;
Average_rate = zeros(length(dLink),number,iter_num+1);
Average_p1 = zeros(length(dLink),number,iter_num+1);
Average_p2 = zeros(length(dLink),number,iter_num+1);
Average_p3 = zeros(length(dLink),number,iter_num+1);
Average_p4 = zeros(length(dLink),number,iter_num+1);

%% Simulation starts
% matlabpool('open');    
for d_L = 1:length( dLink )
dL = dLink( d_L );
fprintf('The %d th calculation starts at dL = %3.2f ... \n',d_L,dL);
[ rate,p1,p2,p3,p4 ] = monte_carlo( TTIs,H,dL,M,iter_num,P1_max,P2_max,P3_max,P4_max,e1,e2,e3,e4,beta1,beta2,beta3,beta4,number,schemes );
Average_rate( d_L,:,:) = rate / TTIs;
Average_p1( d_L,:,:) = p1 / TTIs;
Average_p2( d_L,:,:) = p2 / TTIs;
Average_p3( d_L,:,:) = p3 / TTIs;
Average_p4( d_L,:,:) = p4 / TTIs;
end
% matlabpool('close');

results_folder = './results';
results_file = fullfile(results_folder,...
    sprintf('Chen_GAME_TTI=10_p_max=30_initial=2_sigma=1.mat'));
fprintf('Saving results to ./results/...\n');
save(results_file);
toc;