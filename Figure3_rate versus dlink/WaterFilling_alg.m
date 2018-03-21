function [ pS,pR ] = WaterFilling_alg( H11,H22,H12,H21,beta1,beta2,sigma_n,iter_num,P1_max,P2_max,e1,e2,N )
pS = zeros( 1, iter_num + 1 );
pR = zeros( 1, iter_num + 1 );
muS = zeros( 1, iter_num + 1 );
muR = zeros( 1, iter_num + 1 );
lamdaS = zeros( 1, iter_num + 1 );
lamdaR = zeros( 1, iter_num + 1 );
gammaS = zeros( 1, iter_num + 1 );
gammaR = zeros( 1, iter_num + 1 );

%% Initialize 
a = trace( H12 * H12');
b = trace( H21 * H21');
c = trace( H11 * H11');
d = trace( H22 * H22');
lamda1 = 0.5;
lamda2 = 0.5;
p1 = 1;
p2 = 0.5; 
tau1 = 0.0001;
tau2 = 0.0001;

%% obtain p1£¬p2 with lamda1,lamda2(gamma1,gamma2) fixed
% [ mu1,mu2 ] = get_mu( p1,p2,sigma_n,beta1,beta2,a,b,c,d,P1_max,P2_max );
% gamma1 = lamda1 * ( 1 - beta1 ) * c;
% gamma2 = lamda2 * ( 1 - beta2 ) * d;
% m1 = ( mu1 + gamma1 )^(-1) - ( beta1 * a * p2 + ( sigma_n )^2 ) / ( beta1 * c );   % Initialize p1
% if m1 <= 0
%        p1 = 0;
%     else 
%        if m1 >= P1_max
%           p1 = P1_max;
%        else 
%           p1 = m1;
%        end
% end
% 
% m2 = ( mu2 + gamma2 )^(-1) - ( beta2 * b * p1 + ( sigma_n )^2 ) / ( beta2 * d );   % Initialize p2
% if m2 <= 0
%        p2 = 0;
%     else 
%        if m2 >= P2_max
%           p2 = P2_max;
%        else 
%           p2 = m2;
%        end
% end
pS( :,1 ) = p1;
pR( :,1 ) = p2; 

%% obtain lamda1,lamda2 with p1£¬p2 fixed
% g1 = ( 1 - beta1 ) * c * p1 - e1;
% g2 = ( 1 - beta2 ) * d * p2 - e2;
% n1 = lamda1 - tau1 * g1;              % Initialize lamda1
% n2 = lamda2 - tau2 * g2;              % Initialize lamda2
% if n1 < 0
%      lamda1 = 0;
%  else
%      lamda1 = n1;
% end
% if n2 < 0
%      lamda2 = 0;
%  else
%      lamda2 = n2;
% end
lamdaS ( :,1 ) = lamda1;           
lamdaR ( :,1 ) = lamda2;
  
%% The Interation Loop
for iter = 1 : iter_num
    fprintf( 'The %d th iteration ......',iter );
%% Update p1£¬p2 with lamda1,lamda2 fixed
  fprintf( 'p1 = %3.2f; p2 = %3.2f;...',p1,p2 );
  [ mu1,mu2 ] = get_mu( p1,p2,sigma_n,beta1,beta2,a,b,c,d,P1_max,P2_max );
  muS( :,iter ) = mu1;
  muR( :,iter ) = mu2;
  gammaS( :,iter ) = lamda1 * ( 1 - beta1 ) * c;
  gammaR( :,iter ) = lamda2 * ( 1 - beta2 ) * d;
  R1 = ( muS( :,iter ) + gammaS( :,iter ) )^(-1) - ( beta1 * a * pR( :,iter ) + sigma_n ) / ( beta1 * c );
    if R1 <= 0
       pS( :,iter + 1 ) = 0;
    else 
       if R1 >= P1_max
          pS( :,iter + 1 ) = P1_max;
       else 
          pS( :,iter + 1 ) = R1;
       end
    end
    
   R2 = ( muR( :,iter ) + gammaR( :,iter ) )^(-1) - ( beta2 * b * pS( :,iter +1 ) + sigma_n ) / ( beta2 * d );
    if R2 <= 0
       pR( :,iter + 1 ) = 0;
    else
       if R2 >= P2_max
          pR( :,iter + 1 ) = P2_max;
       else
          pR( :,iter + 1 ) = R2;
       end
    end
    p1 = pS( :,iter + 1 );
    p2 = pR( :,iter + 1 );
%% Update lamda1,lamda2 with p1£¬p2 fixed
   fprintf( 'lamda1 = %5.4f; lamda2 = %5.4f;\n',lamda1,lamda2 );
   g1( :,iter ) = ( 1 - beta1 ) * c * pS( :,iter ) - e1;
   g2( :,iter ) = ( 1 - beta2 ) * d * pR( :,iter ) - e2;
   lamdaS( :,iter + 1 ) = lamdaS( :,iter ) - tau1 * g1( :,iter );
   lamdaR( :,iter + 1 ) = lamdaR( :,iter ) - tau2 * g2( :,iter );
       if lamdaS( :,iter ) - tau1 * g1( :,iter ) < 0
            lamdaS( :,iter + 1 ) = 0;
         else
            lamdaS( :,iter + 1 ) = lamdaS( :,iter ) - tau1 * g1( :,iter );
       end
       if lamdaR( :,iter ) - tau2 * g2( :,iter ) < 0
            lamdaR( :,iter + 1 ) = 0;
         else
            lamdaR( :,iter + 1 ) = lamdaR( :,iter ) - tau2 * g2( :,iter );
       end
   lamda1 = lamdaS( :,iter + 1 );
   lamda2 = lamdaR( :,iter + 1 );
end
end
