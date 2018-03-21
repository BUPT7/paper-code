function [ rate,p1,p2,p3,p4,p]  = monte_carlo( TTIs,H,dL,N,M,iter_num,P1_max,P2_max,P3_max,P4_max,aL,beta1,beta2,beta3,beta4,number,schemes,e )

%%  parameters
rate = zeros( 1,number,iter_num + 1 );
p1 = zeros( 1,number,iter_num + 1 );
p2 = zeros( 1,number,iter_num + 1 );
p3 = zeros( 1,number,iter_num + 1 );
p4 = zeros( 1,number,iter_num + 1 );
p = zeros( 1,number,iter_num + 1 );
sigma_n = 0.01;
d = 1;

%%  Channel Realizations
di1 = sqrt( dL^2 + d^2 );
di2 = sqrt( (2*dL)^2 + d^2 );
di3 = sqrt( (3*dL)^2 + d^2 );
di4 = sqrt( (4*dL)^2 + d^2 );
for ch = 1 : TTIs
    fprintf( 'The %d th TTI...\n',ch );
    H11 = H( :,1:M,ch );
    H22 = H( :,(M+1):(2*M),ch );
    H33 = H( :,(2*M+1):(3*M),ch );
    H44 = H( :,(3*M+1):(4*M),ch );
%     H55 = H( :,(4*M+1):(5*M),ch );
    
    %H12 = sqrt( 1/di1 )*H( :,(2*M+1):(3*M),ch );
    H12 = sqrt( 1/di1 ) * H( :,1:M,ch );
    H13 = sqrt( 1/di2 ) * H( :,1:M,ch );
    H14 = sqrt( 1/di3 ) * H( :,1:M,ch );
%     H15 = sqrt( 1/di4 ) * H( :,1:M,ch );
    
    %H21 = sqrt( 1/di1 )*H( :,(3*M+1):(4*M),ch );
    H21 = sqrt( 1/di1 ) * H( :,(M+1):(2*M),ch );
    H23 = sqrt( 1/di1 ) * H( :,(M+1):(2*M),ch );
    H24 = sqrt( 1/di2 ) * H( :,(M+1):(2*M),ch );
%     H25 = sqrt( 1/di3 ) * H( :,(M+1):(2*M),ch );
    
    H31 = sqrt( 1/di2 ) * H( :,(2*M+1):(3*M),ch );
    H32 = sqrt( 1/di1 ) * H( :,(2*M+1):(3*M),ch );
    H34 = sqrt( 1/di1 ) * H( :,(2*M+1):(3*M),ch );
%     H35 = sqrt( 1/di2 ) * H( :,(2*M+1):(3*M),ch );
    
    H41 = sqrt( 1/di3 ) * H( :,(3*M+1):(4*M),ch );
    H42 = sqrt( 1/di2 ) * H( :,(3*M+1):(4*M),ch );
    H43 = sqrt( 1/di1 ) * H( :,(3*M+1):(4*M),ch );
%     H45 = sqrt( 1/di1 ) * H( :,(3*M+1):(4*M),ch );
    
%     H51 = sqrt( 1/di4 ) * H( :,(4*M+1):(5*M),ch );
%     H52 = sqrt( 1/di3 ) * H( :,(4*M+1):(5*M),ch );
%     H53 = sqrt( 1/di2 ) * H( :,(4*M+1):(5*M),ch );
%     H54 = sqrt( 1/di1 ) * H( :,(4*M+1):(5*M),ch );
    
    
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
        [ p1_id,p2_id,p3_id] = identity( H11,H22,H33,H12,H13,H21,H23,H31,H32,beta1,beta2,beta3,sigma_n,P1_max,P2_max,P3_max,e,aL);
        for iter = 1 : iter_num + 1
            p1_temp = p1_id;
            p2_temp = p2_id;
            p3_temp = p3_id;
           % p4_temp = p4_id;
            p1( 1,1,iter ) = p1( 1,1,iter ) + p1_temp;
            p2( 1,1,iter ) = p2( 1,1,iter ) + p2_temp;
            p3( 1,1,iter ) = p3( 1,1,iter ) + p3_temp;
            %p4( 1,1,iter ) = p4( 1,1,iter ) + p4_temp;
            p( 1,1,iter ) = (p1( 1,1,iter )+p2( 1,1,iter )+ p3( 1,1,iter ))/3;
            [ rate1_temp ] = get_rate2( p1_temp,p2_temp,p3_temp,H11,H22,H33,H12,H13,H21,H23,H31,H32,sigma_n,beta1,beta2,beta3);
            rate( 1,1,iter ) = rate( 1,1,iter ) + rate1_temp;
        end
end
 if schemes(1,2)
      [ p1_wf,p2_wf,p3_wf,p4_wf ] = identity4( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,e,aL );
    for iter = 1 : iter_num + 1
      p1_temp3 = p1_wf;
      p2_temp3 = p2_wf;
      p3_temp3 = p3_wf;
      p4_temp3 = p4_wf;
      p1( 1,2,iter ) = p1( 1,2,iter ) + p1_temp3;
      p2( 1,2,iter ) = p2( 1,2,iter ) + p2_temp3;
      p3( 1,2,iter ) = p3( 1,2,iter ) + p3_temp3;
      p4( 1,2,iter ) = p4( 1,2,iter ) + p4_temp3;
      p( 1,2,iter ) = (p1( 1,2,iter ) + p2( 1,2,iter ) + p3( 1,2,iter ) + p4( 1,2,iter ))/4;
      [ rate2_temp ] = get_rate1( p1_temp3,p2_temp3,p3_temp3,p4_temp3,H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
      rate( 1,2,iter ) = rate( 1,2,iter ) + rate2_temp;
    end
 end

%%  Alogrithm 2, Sequential IWFA
%  if schemes(1,1)
%       [ p1_wf,p2_wf ] = WaterFilling_alg( H11,H22,H12,H21,beta1,beta2,sigma_n,iter_num,P1_max,P2_max,e,aL );
%     for iter = 1 : iter_num + 1
%       p1_temp = p1_wf( :,iter );
%       p2_temp = p2_wf( :,iter );
%       p1( 1,1,iter ) = p1( 1,1,iter ) + p1_temp;
%       p2( 1,1,iter ) = p2( 1,1,iter ) + p2_temp;
%       p( 1,1,iter ) = (p1( 1,1,iter )+p2( 1,1,iter ))/2;
%       [ rate2_temp ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2 );
%       rate( 1,1,iter ) = rate( 1,1,iter ) + rate2_temp;
%     end
%  end
  %  Alogrithm 3, Sequential IWFA,3 users
 if schemes(1,3)
      [ p1_wf1,p2_wf1,p3_wf1 ] = WaterFilling_alg_3(H11,H22,H33,H12,H13,H21,H23,H31,H32,beta1,beta2,beta3,sigma_n,iter_num,P1_max,P2_max,P3_max,e,aL);
    for iter = 1 : iter_num + 1
      p1_temp1 = p1_wf1( :,iter );
      p2_temp1 = p2_wf1( :,iter );
      p3_temp1 = p3_wf1( :,iter );
      
      p1( 1,3,iter ) = p1( 1,3,iter ) + p1_temp1;
      p2( 1,3,iter ) = p2( 1,3,iter ) + p2_temp1;
      p3( 1,3,iter ) = p3( 1,3,iter ) + p3_temp1;
      p( 1,3,iter ) = (p1( 1,3,iter ) + p2( 1,3,iter ) + p3( 1,3,iter ))/3;
      [ rate3_temp ] = get_rate2( p1_temp1,p2_temp1,p3_temp1,H11,H22,H33,H12,H13,H21,H23,H31,H32,sigma_n,beta1,beta2,beta3 );
      rate( 1,3,iter ) = rate( 1,3,iter ) + rate3_temp;
    end
 end
 %%  Alogrithm 4, Sequential IWFA,4 users
 if schemes(1,4)
      [ p1_wf2,p2_wf2,p3_wf2,p4_wf2 ] = WaterFilling_alg_4( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,e,aL );
    for iter = 1 : iter_num + 1
      p1_temp2 = p1_wf2( :,iter );
      p2_temp2 = p2_wf2( :,iter );
      p3_temp2 = p3_wf2( :,iter );
      p4_temp2 = p4_wf2( :,iter );
      p1( 1,4,iter ) = p1( 1,4,iter ) + p1_temp2;
      p2( 1,4,iter ) = p2( 1,4,iter ) + p2_temp2;
      p3( 1,4,iter ) = p3( 1,4,iter ) + p3_temp2;
      p4( 1,4,iter ) = p4( 1,4,iter ) + p4_temp2;
      p( 1,4,iter ) = (p1( 1,4,iter ) + p2( 1,4,iter ) + p3( 1,4,iter ) + p4( 1,4,iter ))/4;
      [ rate4_temp ] = get_rate1( p1_temp2,p2_temp2,p3_temp2,p4_temp2,H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
      rate( 1,4,iter ) = rate( 1,4,iter ) + rate4_temp;
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