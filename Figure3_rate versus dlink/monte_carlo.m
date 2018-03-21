function [ rate,p1,p2]  = monte_carlo( TTIs,H,dL,N,M,iter_num,P1_max,P2_max,e1,e2,beta1,beta2,number,schemes )

%%  parameters
rate = zeros( 1,number,iter_num + 1 );
p1 = zeros( 1,number,iter_num + 1 );
p2 = zeros( 1,number,iter_num + 1 );
sigma_n = 1;
d = 1;

%%  Channel Realizations
di = sqrt( dL^2 + d^2 );
for ch = 1 : TTIs
    fprintf( 'The %d th TTI...\n',ch );
    H11 = H( :,1:M,ch );
    H22 = H( :,(M+1):(2*M),ch );
    H12 = sqrt( 1/di )*H( :,(2*M+1):(3*M),ch );
    H21 = sqrt( 1/di )*H( :,(3*M+1):(4*M),ch );
%     H12 = ( 1/(di^3) )*H( :,(2*M+1):(3*M),ch );
%     H21 = ( 1/(di^3) )*H( :,(3*M+1):(4*M),ch );
%%  Alogrithm 1, Random
%  if schemes(1,1)
%       [ p1_random,p2_random ] = Random( beta1,beta2,H11,H22,e1,e2,iter_num );
%     for iter = 1 : iter_num + 1
%       p1_temp = p1_random( :,iter );
%       p2_temp = p2_random( :,iter );
%       [ rate1_temp ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2 );
%       rate( 1,1,iter ) = rate( 1,1,iter ) + rate1_temp;
%     end 
%  end   
%%  Alogrithm 1, Identity
if schemes(1,1)
        [ p1_id,p2_id ] = identity( H11,H22,H12,H21,beta1,beta2,sigma_n,P1_max,P2_max );
        for iter = 1 : iter_num + 1
            p1_temp = p1_id;
            p2_temp = p2_id;
            [ rate1_temp ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2 );
            rate( 1,1,iter ) = rate( 1,1,iter ) + rate1_temp;
        end
end

%%  Alogrithm 2, Sequential IWFA
 if schemes(1,2)
      [ p1_wf,p2_wf ] = WaterFilling_alg( H11,H22,H12,H21,beta1,beta2,sigma_n,iter_num,P1_max,P2_max,e1,e2,N );
    for iter = 1 : iter_num + 1
      p1_temp = p1_wf( :,iter );
      p2_temp = p2_wf( :,iter );
      p1( 1,2,iter ) = p1( 1,2,iter ) + p1_temp;
      p2( 1,2,iter ) = p2( 1,2,iter ) + p2_temp;
      [ rate2_temp ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2 );
      rate( 1,2,iter ) = rate( 1,2,iter ) + rate2_temp;
    end
 end
%  %%  Alogrithm 3, SDP
%  if schemes(1,3)
%       [ p1_sdp,p2_sdp ] = Filter_SDP( H11,H22,H12,H21,beta1,beta2,sigma_n,iter_num,P1_max,P2_max,e1,e2);
%     for iter = 1 : iter_num + 1
%       p1_temp = p1_sdp( :,iter );
%       p2_temp = p2_sdp( :,iter );
%       [ rate3_temp ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2 );
%       rate( 1,3,iter ) = rate( 1,3,iter ) + rate3_temp;
%     end 
%  end   
end
end