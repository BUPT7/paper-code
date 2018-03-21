% load('./results/Chen_GAME_dL=12_TTI=10_p=20_sigma=1_p1=1_p2=0.5.mat');
% figure(1);
% dLi= 7;
% number = 2;
% init_point = 1;
% line_dot = {'p','-','k';'s','-','k'};
% scheme_name = {'unaided scheme','proposed IWFA scheme'};
% for iter_rate = init_point:number
%     str = strcat( line_dot( iter_rate,1 ),line_dot( iter_rate,2 ),line_dot( iter_rate,3 ) );
%     rate1 = zeros( number - init_point + 1, iter_num + 1 );
%     rate1( :,: ) = Average_rate( dLi,:,: ); 
%     for iter = 1 : iter_num + 1
%         rate1(iter_rate,iter) = Average_rate(dLi,iter_rate,iter);  
%     end
%     semilogy( 0:iter_num,rate1( iter_rate,: ),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3)));
%     hold on;
% end
% xlabel('Iterations');                                        %Label for x-axis
% ylabel('sum rate');                                          %Label for y-axis
% grid on;
% legend(scheme_name(1,1:number),0);
% title('Sum Rate Curve with Iterations');


load('./results/Chen_GAME_dL=14_TTI=10_p=30_sigma=1_p1=1_p2=0.5.mat');
figure(2);
iter = 30;
number = 2;
init_point = 1;
line_dot = {'o',':','k';'s',':','k'};
scheme_name = {'unaided scheme_p=20_sigma=0.5','proposed IWFA scheme_p=20_sigma=0.5'};
for iter_rate = init_point:number
    str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
    rate2 = zeros( length(dLink),number - init_point + 1 );
    for d_L = 1 : length(dLink)
        rate2(d_L,iter_rate) = Average_rate(d_L,iter_rate,iter); 
    end    
    semilogy( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%      plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
    hold on;
end
xlabel('dL');                                                 %Label for x-axis
ylabel('sum rate');                                           %Label for y-axis
grid on;
legend(scheme_name(1,1:number),0);
title('sum rate Curve with dL');
%

% load('./results/Chen_GAME_dL=14_TTI=10_p=20_sigma=0.5_p1=1_p2=0.5.mat');
% number = 2;
% init_point = 1;
% iter = 30;
% line_dot = {'o','-','k';'s','-','k'};
% scheme_name = {'unaided scheme_p=30_sigma=0.5','proposed IWFA scheme_p=30_sigma=0.5'};
% for iter_rate = init_point:number
%     str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
%     rate2 = zeros( length(dLink),number - init_point + 1 );
%     for d_L = 1 : length(dLink)
%         rate2(d_L,iter_rate) = Average_rate(d_L,iter_rate,iter); 
%     end    
%     semilogy( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%     plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );    
% end
% legend(scheme_name(1,1:number),0);
% 
% 
% 
% 
% 
% load('./results/Chen_GAME_dL=14_TTI=10_p=30_sigma=1_p1=1_p2=0.5.mat');
% number = 2;
% init_point = 1;
% iter = 30;
% line_dot = {'o',':','g';'s',':','g'};
% scheme_name = {'unaided scheme_p=20_sigma=1','proposed IWFA scheme_p=20_sigma=1'};
% for iter_rate = init_point:number
%     str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
%     rate2 = zeros( length(dLink),number - init_point + 1 );
%     for d_L = 1 : length(dLink)
%         rate2(d_L,iter_rate) = Average_rate(d_L,iter_rate,iter); 
%     end    
%     semilogy( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%     plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );   
% end
% legend(scheme_name(1,1:number),0);
% 
% 
% 
% load('./results/Chen_GAME_dL=14_TTI=10_p=20_sigma=1_p1=1_p2=0.5.mat');
% number = 2;
% init_point = 1;
% iter = 30;
% line_dot = {'o','-','g';'s','-','g'};
% scheme_name = {'unaided scheme_p=30_sigma=1','proposed IWFA scheme_p=30_sigma=1'};
% for iter_rate = init_point:number
%     str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
%     rate2 = zeros( length(dLink),number - init_point + 1 );
%     for d_L = 1 : length(dLink)
%         rate2(d_L,iter_rate) = Average_rate(d_L,iter_rate,iter); 
%     end    
%     semilogy( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%     plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );  
% end
% legend(scheme_name(1,1:number),0);

% 
% 
% load('./results/Chen_GAME_TTI=10_p=30_sigma=0.1_p1=1_p2=0.5.mat');
% figure(2);
% iter = 30;
% line_dot = {'o','-','r';'s','-','r';};
% for iter_rate = init_point:number
%     str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
%     rate2 = zeros( length(dLink),number - init_point + 1 );
%     for d_L = 1 : length(dLink)
%         rate2(d_L,iter_rate) = Average_rate(d_L,iter_rate,iter); 
%     end    
%     semilogy( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%      plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%     hold on;
% end
% xlabel('dLink');                                              %Label for x-axis
% ylabel('sum rate');                                           %Label for y-axis
% grid on;
% legend(scheme_name(1,1:number),0);
% title('sum rate Curve with dL');
% 
% 
% 
% load('./results/Chen_GAME_TTI=10_p=20_sigma=0.1_p1=1_p2=0.5.mat');
% figure(2);
% iter = 30;
% line_dot = {'o',':','r';'s',':','r';};
% for iter_rate = init_point:number
%     str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
%     rate2 = zeros( length(dLink),number - init_point + 1 );
%     for d_L = 1 : length(dLink)
%         rate2(d_L,iter_rate) = Average_rate(d_L,iter_rate,iter); 
%     end    
%     semilogy( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%      plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%     hold on;
% end
% xlabel('dLink');                                              %Label for x-axis
% ylabel('sum rate');                                           %Label for y-axis
% grid on;
% legend(scheme_name(1,1:number),0);
% title('sum rate Curve with dL');
