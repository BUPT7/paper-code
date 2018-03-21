function [ pS,pR,pO,pQ ] = WaterFilling_alg_6( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,eL )
pS = zeros( 1, iter_num + 1 );
pR = zeros( 1, iter_num + 1 );
pO = zeros( 1, iter_num + 1 );
pQ = zeros( 1, iter_num + 1 );
muS = zeros( 1, iter_num + 1 );
muR = zeros( 1, iter_num + 1 );
muO = zeros( 1, iter_num + 1 );
muQ = zeros( 1, iter_num + 1 );
gammaS = zeros( 1, iter_num + 1 );
gammaR = zeros( 1, iter_num + 1 );
gammaO = zeros( 1, iter_num + 1 );
gammaQ = zeros( 1, iter_num + 1 );
lamdaS = zeros( 1, iter_num + 1 );
lamdaR = zeros( 1, iter_num + 1 );
lamdaO = zeros( 1, iter_num + 1 );
lamdaQ = zeros( 1, iter_num + 1 );

%% Initialize 
p1 = 1;
p2 = 1;
p3 = 1;
p4 = 1; 
pS( :,1 ) = p1;
pR( :,1 ) = p2;
pO( :,1 ) = p3;
pQ( :,1 ) = p4;
lamda1 = 0.7;
lamda2 = 0.7;
lamda3 = 0.7;
lamda4 = 0.7;
lamdaS ( :,1 ) = lamda1;
lamdaR ( :,1 ) = lamda2;
lamdaO ( :,1 ) = lamda3;
lamdaQ ( :,1 ) = lamda4;
tau1 = 0.1;
tau2 = 0.1;
tau3 = 0.1;
tau4 = 0.1;

%% obtain p1£¬p2 with lamda1,lamda2(gamma1,gamma2) fixed
% p1 = 1 \ ( omega1 + gamma1 ) - ( beta1 * norm( H12 * H12' ) * p2 + ( sigma_n )^2 ) \ ( beta1 * norm( H11 * H11' ) );
% p2 = 1 \ ( omega2 + gamma2 ) - ( beta2 * norm( H21 * H21' ) * p1 + ( sigma_n )^2 ) \ ( beta2 * norm( H22 * H22' ) );
% omega1 = ( beta1 * norm( H11 * H11' ) ) \ ( beta1 * norm( H11 * H11') * p1 + beta1 * norm( H12 * H12')  * p2 + ( sigma_n )^2 );
% omega2 = ( beta2 * norm( H22 * H22' ) ) \ ( beta2 * norm( H21 * H21') * p1 + beta2 * norm( H22 * H22')  * p2 + ( sigma_n )^2 );

%% obtain lamda1,lamda2 with p1£¬p2 fixed
% g1 = ( 1 - bata1 ) * norm( H11 * H11' ) * p1 - e1;
% g2 = ( 1 - bata2 ) * norm( H22 * H22' ) * p2 - e2;
% lamda1 = lamda1 - tau1 * g1;
% lamda2 = lamda2 - tau2 * g2;
% gamma1 = lamda1 * ( 1 - beta1 ) * norm( H11 * H11' );
% gamma2 = lamda2 * ( 1 - beta2 ) * norm( H22 * H22' );

%% The Interation Loop
for iter = 1 : iter_num
    fprintf( 'The %d th iteration ......',iter );
    fprintf( 'p1 = %3.2f; p2 = %3.2f;p3 = %3.2f; p4 = %3.2f;...',p1,p2,p3,p4 )
%% Update p1£¬p2 with lamda1,lamda2 fixed
  [ mu1,mu2,mu3,mu4 ] = get_mu4( p1,p2,p3,p4,sigma_n,beta1,beta2,beta3,beta4,H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,P1_max,P2_max,P3_max,P4_max );
  fprintf( 'mu1 = %5.4f; mu2 = %5.4f;mu3 = %5.4f; mu4 = %5.4f;...',mu1,mu2,mu3,mu4 );
  muS( :,iter ) = mu1;
  muR( :,iter ) = mu2;
  muO( :,iter ) = mu3;
  muQ( :,iter ) = mu4;
  gammaS( :,iter ) = lamdaS( :,iter ) * ( 1 - beta1 ) * trace( H11 * H11');
  gammaR( :,iter ) = lamdaR( :,iter ) * ( 1 - beta2 ) * trace( H22 * H22');
  gammaO( :,iter ) = lamdaO( :,iter ) * ( 1 - beta3 ) * trace( H33 * H33');
  gammaQ( :,iter ) = lamdaQ( :,iter ) * ( 1 - beta4 ) * trace( H44 * H44');
  A1 = ( muS( :,iter ) + gammaS( :,iter ) )^(-1) - ( beta1 * trace( H21 * H21') * pR( :,iter ) + beta1 * trace( H31 * H31') * pO( :,iter ) + beta1 * trace( H41 * H41') * pQ( :,iter ) + sigma_n ) / ( beta1 * trace( H11 * H11') );
    if A1 <= 0
       pS( :,iter + 1 ) = 0;
    else 
       if A1 >= P1_max
          pS( :,iter + 1 ) = P1_max;
       else 
          pS( :,iter + 1 ) = A1;
       end
    end
    
    
 A2 = ( muR( :,iter ) + gammaR( :,iter ) )^(-1) - ( beta2 * trace( H12 * H12') * pS( :,iter + 1 ) + beta2 * trace( H32 * H32') * pO( :,iter ) + beta2 * trace( H42 * H42') * pQ( :,iter ) + sigma_n ) / ( beta2 * trace( H22 * H22') );
    if A2 <= 0
       pR( :,iter + 1 ) = 0;
    else
       if A2 >= P2_max
          pR( :,iter + 1 ) = P2_max;
       else
          pR( :,iter + 1 ) = A2;
       end
    end
    
    
  A3 = ( muO( :,iter ) + gammaO( :,iter ) )^(-1) - ( beta3 * trace( H13 * H13') * pS( :,iter + 1 ) + beta3 * trace( H23 * H23') * pR( :,iter + 1 ) + beta3 * trace( H43 * H43') * pQ( :,iter ) + sigma_n ) / ( beta3 * trace( H33 * H33') );
    if A3 <= 0
       pO( :,iter + 1 ) = 0;
    else 
       if A3 >= P3_max
          pO( :,iter + 1 ) = P3_max;
       else 
          pO( :,iter + 1 ) = A3;
       end
    end
    
    
  A4 = ( muQ( :,iter ) + gammaQ( :,iter ) )^(-1) - ( beta4 * trace( H14 * H14') * pS( :,iter + 1 ) + beta4 * trace( H24 * H24') * pR( :,iter + 1 ) + beta4 * trace( H34 * H34') * pO( :,iter + 1 ) + sigma_n ) / ( beta4 * trace( H44 * H44') );
    if A4 <= 0
       pQ( :,iter + 1 ) = 0;
    else
       if A4 >= P4_max
          pQ( :,iter + 1 ) = P4_max;
       else
          pQ( :,iter + 1 ) = A4;
       end
    end
    p1 = pS( :,iter + 1 );
    p2 = pR( :,iter + 1 );
    p3 = pO( :,iter + 1 );
    p4 = pQ( :,iter + 1 );
%% Update lamda1,lamda2 with p1£¬p2 fixed
   fprintf( 'lamda1 = %5.4f; lamda2 = %5.4f;lamda3 = %5.4f; lamda4 = %5.4f;\n',lamda1,lamda2,lamda3,lamda4 );
   g1( :,iter ) = ( 1 - beta1 ) * (trace( H11 * H11') * pS( :,iter ) + trace( H21 * H21') * pR( :,iter ) + trace( H31 * H31') * pO( :,iter ) + trace( H41 * H41') * pQ( :,iter )) - eL;
   g2( :,iter ) = ( 1 - beta2 ) * (trace( H22 * H22') * pR( :,iter)  + trace( H12 * H12') * pS( :,iter ) + trace( H32 * H32') * pO( :,iter ) + trace( H42 * H42') * pQ( :,iter )) - eL;
   g3( :,iter ) = ( 1 - beta3 ) * (trace( H33 * H33') * pO( :,iter ) + trace( H13 * H13') * pS( :,iter ) + trace( H23 * H23') * pR( :,iter ) + trace( H43 * H43') * pQ( :,iter )) - eL;
   g4( :,iter ) = ( 1 - beta4 ) * (trace( H44 * H44') * pO( :,iter ) + trace( H14 * H14') * pS( :,iter ) + trace( H24 * H24') * pR( :,iter ) + trace( H34 * H34') * pO( :,iter )) - eL;
   lamdaS( :,iter + 1 ) = lamdaS( :,iter ) - tau1 * g1( :,iter );
   lamdaR( :,iter + 1 ) = lamdaR( :,iter ) - tau2 * g2( :,iter );
   lamdaO( :,iter + 1 ) = lamdaO( :,iter ) - tau3 * g3( :,iter );
   lamdaQ( :,iter + 1 ) = lamdaQ( :,iter ) - tau4 * g4( :,iter );
   B1 = lamdaS( :,iter ) - tau1 * g1( :,iter );     
       if B1 < 0
            lamdaS( :,iter + 1 ) = 0;
         else
            lamdaS( :,iter + 1 ) = B1;
       end
  
   B2 = lamdaR( :,iter ) - tau2 * g2( :,iter );
       if B2 < 0
            lamdaR( :,iter + 1 ) = 0;
         else
            lamdaR( :,iter + 1 ) = B2;
       end
    
   B3 = lamdaO( :,iter ) - tau3 * g3( :,iter );
       if B3 < 0
            lamdaO( :,iter + 1 ) = 0;
         else
            lamdaO( :,iter + 1 ) = B3;
       end
    
   B4 = lamdaQ( :,iter ) - tau4 * g4( :,iter );
       if B4 < 0
            lamdaQ( :,iter + 1 ) = 0;
         else
            lamdaQ( :,iter + 1 ) = B4;
       end
   lamda1 = lamdaS( :,iter + 1 );
   lamda2 = lamdaR( :,iter + 1 );
   lamda3 = lamdaO( :,iter + 1 );
   lamda4 = lamdaQ( :,iter + 1 );
end
end