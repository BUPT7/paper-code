% load('./results/Chen_GAME_dL=14_TTI=10_p=30_sigma=1_p1=1_p2=0.5.mat'); 
% figure(3);
% eLi = 10;
% number = 2;
% init_point = 1;
% line_dot = {'s','-','b';'o','-','b';};
% scheme_name = {'p1 =1.1','p2 = 0.8'};
% plot_matrix = zeros( 2,iter_num + 1 );
% plot_matrix( 1,: ) = Average_p1( eLi,2,: ); 
% plot_matrix( 2,: ) = Average_p2( eLi,2,: );  
% for iter_p = init_point:number
%     str = strcat( line_dot( iter_p,1 ),line_dot( iter_p,2 ),line_dot( iter_p,3 ) );
%      %{
%     p1 = zeros( 1, iter_num + 1 );
%     p1( :,: ) = Average_p1( dLi,2,: ); 
%     p2 = zeros( 1, iter_num + 1 );
%     p2( :,: ) = Average_p2( dLi,2,: ); 
%     for iter = 1 : iter_num + 1
%         p1(iter_p,iter) = Average_p1(dLi,iter_p,iter); 
%         p2(iter_p,iter) = Average_p2(dLi,iter_p,iter);
%     end
%     %}
%     plot( 0:iter_num,plot_matrix( iter_p,: ),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_p,3)));
%     hold on;
% end
% xlabel('Iterations');                                                 %Label for x-axis
% ylabel('power allocations');                                          %Label for y-axis
% grid on;
% legend(scheme_name(1,1:number),0);
% title('power allocations Curve with Iterations');

load('./results/Chen_GAME_dL=14_TTI=10_p=30_sigma=1_p1=1_p2=0.5.mat');
figure(1);
iter = 10;
number = 6;
init_point = 1;
line_dot = {'o',':','r';'s',':','k';'o','-','r';'o','-','k';'o',':','b',;'o','-','b'};
scheme_name = {'unaided scheme (3 links, \alpha = 0.5)','unaided scheme (4 links, \alpha = 0.5)','proposed IWFA scheme (3 links, \alpha = 0.5)','proposed IWFA scheme (4 links, \alpha = 0.5)','proposed IWFA scheme (4 links, \alpha = 0.6)','proposed IWFA scheme (4 links, \alpha = 0.7)'};
for iter_rate = init_point:number
    str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
    p_1 = zeros( length(eLink),number - init_point + 1 );
    p_2 = zeros( length(eLink),number - init_point + 1 );
    for e_L = 1 : length(eLink)
        p(e_L,iter_rate) = Average_p(e_L,iter_rate,iter+1);
        %rate(e_L,iter_rate) = Average_rate(e_L,iter_rate,iter);
        %p_sum(e_L,iter_rate)= p_1(e_L,iter_rate)+p_2(e_L,iter_rate);
    end    
    semilogy( eLink,p(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%      plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
    hold on;
end
xlabel('the harvest energy threshold');                                                 %Label for x-axis
ylabel('Average transmit power');   
%axis([0.1 1 0.25 0.5]);
grid on;
legend(scheme_name(1,1:number),0);
title('Average transmit powe with the harvest energy threshold');

% load('./results/Chen_GAME_POWERallo1_p=30_d=5_0.7_1.3.mat'); 
% figure(5);
% dLi= 7;
% number = 2;
% init_point = 1;
% line_dot = {'s','-','g';'o','-','g'};
% scheme_name = {'p1 = 2','p2 = 2'};
% plot_matrix = zeros( 2,iter_num + 1 );
% plot_matrix( 1,: ) = Average_p1( dLi,2,: ); 
% plot_matrix( 2,: ) = Average_p2( dLi,2,: ); 
% for iter_p = init_point:number
%     str = strcat( line_dot( iter_p,1 ),line_dot( iter_p,2 ),line_dot( iter_p,3 ) );
% %      %{
% %     p1 = zeros( 1, iter_num + 1 );
% %     p1( :,: ) = Average_p1( dLi,2,: ); 
% %     p2 = zeros( 1, iter_num + 1 );
% %     p2( :,: ) = Average_p2( dLi,2,: ); 
% %     for iter = 1 : iter_num + 1
% %         p1(iter_p,iter) = Average_p1(dLi,iter_p,iter); 
% %         p2(iter_p,iter) = Average_p2(dLi,iter_p,iter);
% %     end
% %     %}
%     plot( 0:iter_num,plot_matrix( iter_p,: ),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_p,3)));
%     hold on;
% end
% xlabel('Iterations');                                                 %Label for x-axis
% ylabel('power allocations');                                          %Label for y-axis
% grid on;
% legend(scheme_name(1,1:number),0);
% title('power allocations Curve with Iterations');