function [ rate,p1,p2,p3,p4,p]  = monte_carlo( TTIs,H,dL,N,M,iter_num,P1_max,P2_max,P3_max,P4_max,eL,beta1,beta2,beta3,beta4,number,schemes  )

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
    H11 = H( :,1:M,ch );
    H22 = H( :,(M+1):(2*M),ch );
    H33 = H( :,(2*M+1):(3*M),ch );
    H44 = H( :,(3*M+1):(4*M),ch );

 
%%  Alogrithm 1, Identity
if schemes(1,1)
        [ p1_id,p2_id,p3_id ] = identity( H11,H22,H33,H12,H13,H21,H23,H31,H32,beta1,beta2,beta3,sigma_n,iter_num,P1_max,P2_max,P3_max,eL);  
        for iter = 1 : iter_num + 1
            p1_temp = p1_id;
            p2_temp = p2_id;
            p3_temp = p3_id;
            %p4_temp = p4_id;
            p1( 1,1,iter ) = p1( 1,1,iter ) + p1_temp;
            p2( 1,1,iter ) = p2( 1,1,iter ) + p2_temp;
            p3( 1,1,iter ) = p3( 1,1,iter ) + p3_temp;
            %p4( 1,1,iter ) = p4( 1,1,iter ) + p4_temp;
            p( 1,1,iter ) = (p1( 1,1,iter )+p2( 1,1,iter )+p3( 1,1,iter ))/3;
            [ rate1_temp ] = get_rate2( p1( 1,1,iter ),p2( 1,1,iter ),p3( 1,1,iter ),H11,H22,H33,H12,H13,H21,H23,H31,H32,sigma_n,beta1,beta2,beta3 );
            rate( 1,1,iter ) = rate( 1,1,iter ) + rate1_temp;
         end
end
if schemes(1,2)
        [ p1_id2,p2_id2,p3_id2,p4_id2 ] = identity1(H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,eL);
        for iter = 1 : iter_num + 1
            p1_temp2 = p1_id2;
            p2_temp2 = p2_id2;
            p3_temp2 = p3_id2;
            p4_temp2 = p4_id2;
            p1( 1,2,iter ) = p1( 1,2,iter ) + p1_temp2;
            p2( 1,2,iter ) = p2( 1,2,iter ) + p2_temp2;
            p3( 1,2,iter ) = p3( 1,2,iter ) + p3_temp2;
            p4( 1,2,iter ) = p4( 1,2,iter ) + p4_temp2;
            %p4( 1,1,iter ) = p4( 1,1,iter ) + p4_temp;
            p( 1,2,iter ) = (p2( 1,2,iter )+p2( 1,2,iter )+p3( 1,2,iter )+p4( 1,2,iter ))/4;
            [ rate2_temp ] = get_rate1( p1( 1,2,iter ), p2( 1,2,iter ),p3( 1,2,iter ),p4( 1,2,iter ),H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
            rate( 1,2,iter ) = rate( 1,2,iter ) + rate2_temp;
        end
end

%  Alogrithm 2, Sequential IWFA
%  if schemes(1,2)
%       [ p1_wf,p2_wf,p3_wf,p4_wf ] = WaterFilling_alg(H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,eL);
%     for iter = 1 : iter_num + 1
%       p1_temp2 = p1_wf( :,iter );
%       p2_temp2 = p2_wf( :,iter );
%        p3_temp2 = p3_wf( :,iter );
%         p4_temp2 = p4_wf( :,iter );
%       p1( 1,2,iter ) = p1( 1,2,iter ) + p1_temp2;
%       p2( 1,2,iter ) = p2( 1,2,iter ) + p2_temp2;
%       p3( 1,2,iter ) = p3( 1,2,iter ) + p3_temp2;
%       p4( 1,2,iter ) = p4( 1,2,iter ) + p4_temp2;
%       p( 1,2,iter ) = (p1( 1,2,iter )+p2( 1,2,iter )+p3( 1,2,iter )+p4( 1,2,iter ))/4;
%       [ rate2_temp ] = get_rate1( p1( 1,2,iter ),p2( 1,2,iter ),p3( 1,2,iter ),p4( 1,2,iter ),H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
%       rate( 1,2,iter ) = rate( 1,2,iter ) + rate2_temp;
%     end
%  end
  %%  Alogrithm 3, Sequential IWFA,3 users
 if schemes(1,3)
      [ p1_wf1,p2_wf1,p3_wf1 ] = WaterFilling_alg_3(H11,H22,H33,H12,H13,H21,H23,H31,H32,beta1,beta2,beta3,sigma_n,iter_num,P1_max,P2_max,P3_max,eL);
     for iter = 1 : iter_num + 1
      p1_temp3 = p1_wf1( :,iter );
      p2_temp3 = p2_wf1( :,iter );
      p3_temp3 = p3_wf1( :,iter );
      
      p1( 1,3,iter ) = p1( 1,3,iter ) + p1_temp3;
      p2( 1,3,iter ) = p2( 1,3,iter ) + p2_temp3;
      p3( 1,3,iter ) = p3( 1,3,iter ) + p3_temp3;
      p( 1,3,iter ) = (p1( 1,3,iter ) + p2( 1,3,iter ) + p3( 1,3,iter ))/3;
      p( 1,3,iter ) = ( p1( 1,3,iter ) + p2( 1,3,iter ) + p3( 1,3,iter ))/3;
      [ rate3_temp ] = get_rate2( p1( 1,3,iter ),p2( 1,3,iter ),p3( 1,3,iter ),H11,H22,H33,H12,H13,H21,H23,H31,H32,sigma_n,beta1,beta2,beta3 );
     %[ rate3_temp ] = get_rate2( p1_temp3,p2_temp3,p3_temp3,H11,H22,H33,H12,H13,H21,H23,H31,H32,sigma_n,beta1,beta2,beta3 );
      rate( 1,3,iter ) = rate( 1,3,iter ) + rate3_temp;
     end
 end
 %%  Alogrithm 4, Sequential IWFA,4 users
%  if schemes(1,4)
%       [ p1_wf2,p2_wf2,p3_wf2,p4_wf2 ] = WaterFilling_alg_4( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,eL );
%     for iter = 1 : iter_num + 1
%       p1_temp4 = p1_wf2( :,iter );
%       p2_temp4 = p2_wf2( :,iter );
%       p3_temp4 = p3_wf2( :,iter );
%       p4_temp4 = p4_wf2( :,iter );
%       p1( 1,4,iter ) = p1( 1,4,iter ) + p1_temp4;
%       p2( 1,4,iter ) = p2( 1,4,iter ) + p2_temp4;
%       p3( 1,4,iter ) = p3( 1,4,iter ) + p3_temp4;
%       p4( 1,4,iter ) = p4( 1,4,iter ) + p4_temp4;
%       p( 1,4,iter ) = (p1( 1,4,iter ) + p2( 1,4,iter ) + p3( 1,4,iter ) + p4( 1,4,iter ))/4;
%       [ rate4_temp ] = get_rate1( p1( 1,4,iter ),p2( 1,4,iter ),p3( 1,4,iter ),p4( 1,4,iter ),H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
%       rate( 1,4,iter ) = rate( 1,4,iter ) + rate4_temp;
%     end
%  end
  if schemes(1,4)
      [ p1_wf3,p2_wf3,p3_wf3,p4_wf3 ] = WaterFilling_alg_5( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,eL );
    for iter = 1 : iter_num + 1
      p1_temp5 = p1_wf3( :,iter );
      p2_temp5 = p2_wf3( :,iter );
      p3_temp5 = p3_wf3( :,iter );
      p4_temp5 = p4_wf3( :,iter );
      p1( 1,4,iter ) = p1( 1,4,iter ) + p1_temp5;
      p2( 1,4,iter ) = p2( 1,4,iter ) + p2_temp5;
      p3( 1,4,iter ) = p3( 1,4,iter ) + p3_temp5;
      p4( 1,4,iter ) = p4( 1,4,iter ) + p4_temp5;
      p( 1,4,iter ) = (p1( 1,4,iter ) + p2( 1,4,iter ) + p3( 1,4,iter ) + p4( 1,4,iter ))/4;
      [ rate5_temp ] = get_rate1( p1( 1,4,iter ),p2( 1,4,iter ),p3( 1,4,iter ),p4( 1,4,iter ),H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
      rate( 1,4,iter ) = rate( 1,4,iter ) + rate5_temp;
    end
  end
  if schemes(1,5)
      [ p1_wf4,p2_wf4,p3_wf4,p4_wf4 ] = WaterFilling_alg_6( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,eL );
    for iter = 1 : iter_num + 1
      p1_temp4 = p1_wf4( :,iter );
      p2_temp4 = p2_wf4( :,iter );
      p3_temp4 = p3_wf4( :,iter );
      p4_temp4 = p4_wf4( :,iter );
      p1( 1,5,iter ) = p1( 1,5,iter ) + p1_temp4;
      p2( 1,5,iter ) = p2( 1,5,iter ) + p2_temp4;
      p3( 1,5,iter ) = p3( 1,5,iter ) + p3_temp4;
      p4( 1,5,iter ) = p4( 1,5,iter ) + p4_temp4;
      p( 1,5,iter ) = (p1( 1,5,iter ) + p2( 1,5,iter ) + p3( 1,5,iter ) + p4( 1,5,iter ))/4;
      [ rate6_temp ] = get_rate1( p1( 1,5,iter ),p2( 1,5,iter ),p3( 1,5,iter ),p4( 1,5,iter ),H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
      rate( 1,5,iter ) = rate( 1,5,iter ) + rate6_temp;
    end
  end
   if schemes(1,6)
      [ p1_wf5,p2_wf5,p3_wf5,p4_wf5 ] = WaterFilling_alg_7( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,eL );
    for iter = 1 : iter_num + 1
      p1_temp6 = p1_wf5( :,iter );
      p2_temp6 = p2_wf5( :,iter );
      p3_temp6 = p3_wf5( :,iter );
      p4_temp6 = p4_wf5( :,iter );
      p1( 1,6,iter ) = p1( 1,6,iter ) + p1_temp6;
      p2( 1,6,iter ) = p2( 1,6,iter ) + p2_temp6;
      p3( 1,6,iter ) = p3( 1,6,iter ) + p3_temp6;
      p4( 1,6,iter ) = p4( 1,6,iter ) + p4_temp6;
      p( 1,6,iter ) = (p1( 1,6,iter ) + p2( 1,6,iter ) + p3( 1,6,iter ) + p4( 1,6,iter ))/4;
      [ rate4_temp ] = get_rate1( p1( 1,6,iter ),p2( 1,6,iter ),p3( 1,6,iter ),p4( 1,6,iter ),H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,sigma_n,beta1,beta2,beta3,beta4 );
      rate( 1,6,iter ) = rate( 1,6,iter ) + rate4_temp;
    end
   end
end
end