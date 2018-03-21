function [ rate ] = get_rate2( p1_temp,p2_temp,p3_temp,H11,H22,H33,H12,H13,H21,H23,H31,H32,sigma_n,beta1,beta2,beta3 )

Z1 = ( beta1 * trace( H11 * H11' )* p1_temp ) / ( beta1 * trace( H21 * H21') * p2_temp + beta1 * trace( H31 * H31') * p3_temp  + sigma_n );
rate1 = log2( 1 + Z1 );
Z2 = ( beta2 * trace( H22 * H22') * p2_temp ) / ( beta2 * trace( H12 * H12') * p1_temp + beta2 * trace( H32 * H32') * p3_temp + sigma_n );
rate2 = log2( 1 + Z2 );
Z3 = ( beta3 * trace( H33 * H33' )* p3_temp ) / ( beta3 * trace( H13 * H13') * p1_temp + beta3 * trace( H23 * H23') * p2_temp + sigma_n );
rate3 = log2( 1 + Z3 );

rate = rate1 + rate2 + rate3 ;
end