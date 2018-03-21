clear all;
% absolute distance
D_s1_s2 = 1;
R_s1 = 1;
% ralative distance
d = D_s1_s2 / R_s1;  
channels_folder = './channels';
for k_ = 1:2                       % loop the max Node D's anteenas number M 
    for m_ = 1:3                   % loop TTIs = 10,100,1000
        ttis = 10^m_;
        string1 = strcat('H_1_',num2str(k_),'_',num2str(ttis),'.mat');
        H = init_network(ttis,d,1,k_);
        channel_file = fullfile(channels_folder,sprintf(string1));
        save(channel_file);
    end
end
for k_ = 1:3                       % loop the max Node S,D's anteenas number M 
    for m_ = 1:3                   % loop TTIs = 10,100,1000
        ttis = 10^m_;
        string1 = strcat('H_2_',num2str(k_),'_',num2str(ttis),'.mat');
        H = init_network(ttis,d,2,k_);
        channel_file = fullfile(channels_folder,sprintf(string1));
        save(channel_file);
    end
end
for k_ = 1:4                       % loop the max Node S,D's anteenas number M 
    for m_ = 1:3                   % loop TTIs = 10,100,1000
        ttis = 10^m_;
        string2 = strcat('H_3_',num2str(k_),'_',num2str(ttis),'.mat');
        H = init_network(ttis,d,3,k_);
        channel_file = fullfile(channels_folder,sprintf(string2));
        save(channel_file);
    end
end
