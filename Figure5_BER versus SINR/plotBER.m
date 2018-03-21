load('./results/Chen_GAME_dL=14_TTI=10_p=30_sigma=1_p1=1_p2=0.5.mat');
figure(1);
iter = 10;
number = 2;
init_point = 1;
line_dot = {'o',':','b';'s',':','k';'o','-','r';'o','-','r'};
scheme_name = {'proposed IWFA scheme, 3 users'};
for iter_rate = init_point:number
    str = strcat(line_dot(iter_rate,1),line_dot(iter_rate,2),line_dot(iter_rate,3));
%     p_1 = zeros( length(aLink),number - init_point + 1 );
%     p_2 = zeros( length(aLink),number - init_point + 1 );
    for Eb_No = 1 : length(EbNo)
        errs(Eb_No,iter_rate) = Average_BER(Eb_No,iter_rate,iter);
        %rate(e_L,iter_rate) = Average_rate(e_L,iter_rate,iter);
        %p_sum(e_L,iter_rate)= p_1(e_L,iter_rate)+p_2(e_L,iter_rate);
    end    
    semilogy( EbNo,errs(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
%      plot( dLink,rate2(:,iter_rate),cell2mat(str),'MarkerFaceColor',cell2mat(line_dot(iter_rate,3) ) );
    hold on;
end
xlabel('the pricing factor');                                                 %Label for x-axis
ylabel('sum-rate');   
%axis([0 10 2 6]); %Label for y-axis
grid on;
legend(scheme_name(1,1:number),0);
title('sum rate Curve with the pricing factor');
