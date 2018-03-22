function [ rate ,p1,p2,p3,p4 ]  = monte_carlo( TTIs,H,dL,M,iter_num,P1_max,P2_max,P3_max,P4_max,e1,e2,e3,e4,beta1,beta2,beta3,beta4,number,schemes )

%%  parameters
rate = zeros( 1,number,iter_num + 1 );
p1 = zeros( 1,number,iter_num + 1 );
p2 = zeros( 1,number,iter_num + 1 );
p3 = zeros( 1,number,iter_num + 1 );
p4 = zeros( 1,number,iter_num + 1 );
sigma_n = 1;   %%( sigma_n )^2 
d = 1;

%%  Channel Realizations
di1 = sqrt( dL^2 + d^2 );
di2 = sqrt( (2*dL)^2 + d^2 );
di3 = sqrt( (3*dL)^2 + d^2 );
for ch = 1 : TTIs
    fprintf( 'The %d th TTI...\n',ch );
    H11 = H( :,1:M,ch );
    H22 = H( :,(M+1):(2*M),ch );
    H33 = H( :,(2*M+1):(3*M),ch );
    H44 = H( :,(3*M+1):(4*M),ch );
    H12 = sqrt( 1/di1 ) * H( :,1:M,ch );
    H13 = sqrt( 1/di2 ) * H( :,1:M,ch );
    H14 = sqrt( 1/di3 ) * H( :,1:M,ch );
    H21 = sqrt( 1/di1 ) * H( :,(M+1):(2*M),ch );
    H23 = sqrt( 1/di1 ) * H( :,(M+1):(2*M),ch );   
    H24 = sqrt( 1/di2 ) * H( :,(M+1):(2*M),ch );
    H31 = sqrt( 1/di2 ) * H( :,(2*M+1):(3*M),ch );
    H32 = sqrt( 1/di1 ) * H( :,(2*M+1):(3*M),ch );
    H34 = sqrt( 1/di1 ) * H( :,(2*M+1):(3*M),ch );
    H41 = sqrt( 1/di3 ) * H( :,(3*M+1):(4*M),ch );
    H42 = sqrt( 1/di2 ) * H( :,(3*M+1):(4*M),ch );
    H43 = sqrt( 1/di1 ) * H( :,(3*M+1):(4*M),ch );

% %%  Alogrithm 1, Random
%  if schemes(1,1)
%       [ p1_random,p2_random ] = Random( beta1,beta2,H11,H22,e1,e2,P1_max,P2_max,iter_num );
%     for iter = 1 : iter_num + 1
%       p1_temp = p1_random( :,iter );
%       p2_temp = p2_random( :,iter );
%       [ rate1_temp ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2 );
%       rate( 1,1,iter ) = rate( 1,1,iter ) + rate1_temp;
%     end 
%  end   
 
%%  Alogrithm 1, Identity
if schemes(1,1)
        [ p1_id,p2_id,p3_id,p4_id ] = identity( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,e1,e2,e3,e4 );
        for iter = 1 : iter_num + 1
            p1_temp = p1_id;
            p2_temp = p2_id;
            p3_temp = p3_id;
            p4_temp = p4_id;
            [ rate1_temp ] = get_rate( p1_temp,p2_temp,p3_temp,p4_temp,H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
            rate( 1,1,iter ) = rate( 1,1,iter ) + rate1_temp;
        end
end

%%  Alogrithm 2, Sequential IWFA
 if schemes(1,2)
      [ p1_wf,p2_wf,p3_wf,p4_wf ] = WaterFilling_alg( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,e1,e2,e3,e4 );
    for iter = 1 : iter_num + 1
      p1_temp = p1_wf( :,iter );
      p2_temp = p2_wf( :,iter );
      p3_temp = p3_wf( :,iter );
      p4_temp = p4_wf( :,iter );
      p1( 1,2,iter ) = p1( 1,2,iter ) + p1_temp;
      p2( 1,2,iter ) = p2( 1,2,iter ) + p2_temp;
      p3( 1,2,iter ) = p3( 1,2,iter ) + p3_temp;
      p4( 1,2,iter ) = p4( 1,2,iter ) + p4_temp;
      [ rate2_temp ] = get_rate( p1_temp,p2_temp,p3_temp,p4_temp,H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
      rate( 1,2,iter ) = rate( 1,2,iter ) + rate2_temp;
    end
 end
 %%  Alogrithm 3, Optimal
%  if schemes(1,3)
%       [ p1_opt,p2_opt ] = Optimal_alg( beta1,beta2,H11,H22,e1,e2,P1_max,P2_max,iter_num );
%     for iter = 1 : iter_num + 1
%       p1_temp = p1_opt( :,iter );
%       p2_temp = p2_opt( :,iter );
%       [ rate3_temp ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2 );
%       rate( 1,3,iter ) = rate( 1,3,iter ) + rate3_temp;
%     end 
%  end   
end
end