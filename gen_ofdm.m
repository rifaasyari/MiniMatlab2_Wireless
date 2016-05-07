function [ X_MMSE ] = gen_ofdm(h, SNR,modOrd,mu,N_pts,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%modOrd = 2;             % constellation size = 2^modOrd
EbNoVec = SNR - 10*log10(modOrd);        % Eb/No in dB
%a local random stream to be used by random number generators for
% repeatability.
%hStr = RandStream('mt19937ar');

% Create PSK modulator and demodulator System objects
%hMod   = comm.PSKModulator(...
%            'ModulationOrder',  2^modOrd, ...
%            'PhaseOffset',      0, ...
%            'BitInput',         true);
%hDemod = comm.PSKDemodulator( ...
%            'ModulationOrder',  2^modOrd, ...
%            'PhaseOffset',      0, ...
%            'BitOutput',        true);

        
% M-ary QAM    %qammod(dataSymbolsIn,M,0);     and qamdemod(receivedSignal,M);     
%modOrd                      % Size of signal constellation
k = log2(modOrd);                % Number of bits per symbol
n = k*64*100;                  % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor      
rng default                 % Use default random number generator
dataIn = randi([0 1],n,1);  % Generate vector of binary data
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers

% Pre-allocate variables to store BER results for speed
BER_ZF =  zeros(length(EbNoVec), (1/N)*length(dataSymbolsIn)) ;
BER_PRECODING =  zeros(length(EbNoVec),(1/N)*length(dataSymbolsIn)) ;
BER_MMSE = zeros(length(EbNoVec), (1/N)*length(dataSymbolsIn)) ;
BER_BASELINE = zeros(length(EbNoVec), (1/N)*length(dataSymbolsIn)) ;

dataSymbolsIn_mat = reshape(dataSymbolsIn,[N, 1, (1/N)*length(dataSymbolsIn)]);



%%

H_freq = transpose(fft(h,N_pts)) ; 


H = zeros( N_pts, mu+N_pts) ;
for kk = 1 : N_pts
    H(kk,kk:(kk+mu)) = [h ,zeros(1,mu+1-length(h))] ;
end




% Loop over N symbol messages
for idx = 1:length(SNR) 
    
    for p = 1:1:(1/N)*length(dataSymbolsIn)  
   

    % Calculate SNR from EbNo for each independent transmission link
    snrIndB = SNR(idx);
    snrLinear = 10^(0.1*snrIndB);


    % Create random bit vector to modulate
    msg =  dataSymbolsIn_mat(:,:,p) ; 

    % Modulate data
    txSig = qammod(msg,modOrd,0);

    x = ifft(txSig,N_pts); 
    
    x_prefix =   x(end:-1:end-(mu-1) );
    
    x_input =  [ flipud(x) ;  x_prefix] ; 
    
    %
   
    sig_power = var(x_input,1);      
    % Add noise to faded data   y = H* x + n
    %n = sqrt(sig_power*10^(-snrIndB/10))*(randn(N_pts, 1) + j*randn(N_pts, 1))/sqrt(2) ; 
    %y_n  = H*x_input +  n ;
     output = H*x_input;
    
    y_n = awgn(output,snrIndB, sig_power) ;
    %y_remove_prefix = y_n(mu+1:end) ;
   
    y_ordered = flipud(y_n);
    Y_freq = fft(y_ordered,N_pts) ;
    
    
    
    %BASELINE (NO EQUALIZATION)N
    rx_BASELINE = qamdemod(Y_freq,modOrd); 
    BER_BASELINE(idx,p)  = biterr(rx_BASELINE,msg);
    
    
      
    %ZERO-FORCING EQUALIZATION 
 
    %Demodulate data
    X_ZF = Y_freq ./ H_freq ;
    rxSig_ZF = qamdemod(X_ZF,modOrd);
    BER_ZF(idx,p)= biterr(rxSig_ZF,msg);
    
       
    %MMSE  EQUALIZATION
    N_0 = sig_power*10^(-snrIndB/10); 
    X_MMSE = Y_freq ./ (H_freq + N_0*ones(size(H_freq)));
    rxSig_MMSE = qamdemod(X_MMSE,modOrd);
    BER_MMSE(idx,p) = biterr(rxSig_MMSE,msg) ;  
       
    
     
     
    end
end    





end

