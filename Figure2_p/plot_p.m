load('./results/Chen_GAME_TTI=100_p_max=25_initial=1.1_1.5_1.2_1.3_sigma=1.mat'); 
figure(1);
dLi= 6;
number = 4;
init_point = 1;
line_dot = {'s','-','k';'d','-','k';'^','-','k';'o','-','k'};
scheme_name = {'p1 = 1','p2 = 1','p3 = 1','p4 = 1'};
plot_matrix = zeros( 4,iter_num + 1 );
plot_matrix( 1,: ) = Average_p1( dLi,2,: ); 
plot_matrix( 2,: ) = Average_p2( dLi,2,: ); 
plot_matrix( 3,: ) = Average_p3( dLi,2,: ); 
plot_matrix( 4,: ) = Average_p4( dLi,2,: ); 
for iter_p = init_point:number
    str = strcat( line_dot( iter_p,1 ),line_dot( iter_p,2 ),line_dot( iter_p,3 ) );
    plot( 0:iter_num,plot_matrix( iter_p,: ),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_p,3)));
    hold on;
end
xlabel('Iterations');                                                         %Label for x-axis
ylabel('Power Allocation Strategy');                                          %Label for y-axis
grid on;
% legend(scheme_name(1,1:number),0);
title('Power Allocations Strategy Curve with Iterations');


load('./results/Chen_GAME_TTI=100_p_max=25_initial=1.4_1.3_0.8_1.1_sigma=1.mat'); 
figure(1);
dLi= 6;
number = 4;
init_point = 1;
line_dot = {'*',':','k';'.',':','k';'x',':','k';'+',':','k'};
scheme_name = {'p1 = 2','p2 = 2','p3 = 2','p4 = 2'};
plot_matrix = zeros( 4,iter_num + 1 );
plot_matrix( 1,: ) = Average_p1( dLi,2,: ); 
plot_matrix( 2,: ) = Average_p2( dLi,2,: ); 
plot_matrix( 3,: ) = Average_p3( dLi,2,: ); 
plot_matrix( 4,: ) = Average_p4( dLi,2,: );  
for iter_p = init_point:number
    str = strcat( line_dot( iter_p,1 ),line_dot( iter_p,2 ),line_dot( iter_p,3 ) );
    plot( 0:iter_num,plot_matrix( iter_p,: ),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_p,3)));
    hold on;
end
xlabel('Iterations');                                                 %Label for x-axis
ylabel('Power Allocation Strategy');                                          %Label for y-axis
grid on;
legend(scheme_name(1,1:number),0);
title('Power Allocation Strategy Curve with Iterations');