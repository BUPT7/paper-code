function [ rate ] = get_rate( p1_temp,p2_temp,H11,H22,H12,H21,sigma_n,beta1,beta2  )

Z1 = ( beta1 * trace( H11 * H11' )* p1_temp ) / ( beta1 * trace( H21 * H21') * p2_temp + sigma_n );
rate1 = log2( 1 + Z1 );
Z2 = ( beta2 * trace( H22 * H22') * p2_temp ) / ( beta2 * trace( H12 * H12') * p1_temp + sigma_n );
rate2 = log2( 1 + Z2 );
rate = rate1 + rate2;
end

