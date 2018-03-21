function [ H ] = init_network( TTIs,d,N,M )
H = zeros(N,4*M,length(TTIs));

    for T = 1:1:TTIs
        H1 = sqrt(1/d)*(randn(N,M)+1i*randn(N,M))/sqrt(2);% the channel of S1 to D1
        H2 = sqrt(1/d)*(randn(N,M)+1i*randn(N,M))/sqrt(2);% the channel of S2 to D2
%         H1 = (1/(d^3))*(randn(N,M)+1i*randn(N,M))/sqrt(2);% the channel of S1 to D1
%         H2 = (1/(d^3))*(randn(N,M)+1i*randn(N,M))/sqrt(2);% the channel of S2 to D2
        H3 = (randn(N,M)+1i*randn(N,M))/sqrt(2);% the iterference channel of S1 to D2
        H4 = (randn(N,M)+1i*randn(N,M))/sqrt(2);% the iterference channel of S2 to D1
        H(:,:,T) = [ H1 H2 H3 H4 ];    
    end
end



