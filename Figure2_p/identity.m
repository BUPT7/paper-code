function [ p1,p2,p3,p4 ] = idetity( H11,H22,H33,H44,H12,H13,H14,H21,H23,H24,H31,H32,H34,H41,H42,H43,beta1,beta2,beta3,beta4,sigma_n,iter_num,P1_max,P2_max,P3_max,P4_max,e1,e2,e3,e4 );
  p1 = 1.0;
  p2 = 1.6;
  p3 = 0.8;
  p4 = 0.5;
%   lamda1 = 0.5;
%   lamda2 = 0.5;
%   a = trace( H12 * H12');
%   b = trace( H21 * H21');
%   c = trace( H11 * H11');
%   d = trace( H22 * H22');
%   [ mu1,mu2 ] = get_mu( p1,p2,sigma_n,beta1,beta2,a,b,c,d,P1_max,P2_max );
%   gamma1 = lamda1 * ( 1 - beta1 ) * c;
%   gamma2 = lamda2 * ( 1 - beta2 ) * d;
%   m1 = ( mu1 + gamma1 )^(-1) - ( beta1 * a * p2 + ( sigma_n )^2 ) / ( beta1 * c );
%     if m1 <= 0
%        p1 = 0;
%     else 
%        if m1 >= P1_max
%           p1 = P1_max;
%        else 
%           p1 = m1;
%        end
%     end
%     
%    m2 = ( mu2 + gamma2 )^(-1) - ( beta2 * b * p1 + ( sigma_n )^2 ) / ( beta2 * d );
%     if m2 <= 0
%        p2 = 0;
%     else
%        if m2 >= P2_max
%           p2 = P2_max;
%        else
%           p2 = m2;
%        end
%     end
    fprintf( 'p1 = %3.2f; p2 = %3.2f;p3 = %3.2f; p4 = %3.2f;...\n',p1,p2,p3,p4);
end

