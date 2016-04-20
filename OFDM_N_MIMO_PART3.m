function [ mu_ZF , mu_MMSE , mu_NO_EQ ] = OFDM_N_MIMO_PART3(h, SNR,modOrd,mu,M_r,M_t,N_pts,N)
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
BER_ZF =  zeros(length(EbNoVec), (1/N)*length(dataSymbolsIn) , N_pts) ;
BER_PRECODING =  zeros(length(EbNoVec),(1/N)*length(dataSymbolsIn), N_pts) ;
BER_MMSE = zeros(length(EbNoVec), (1/N)*length(dataSymbolsIn), N_pts) ;
BER_BASELINE = zeros(length(EbNoVec), (1/N)*length(dataSymbolsIn),N_pts) ;

dataSymbolsIn_mat = reshape(dataSymbolsIn,[N, 1, (1/N)*length(dataSymbolsIn)]);



%%

H_freq_11 = (fft(h(1,:),N_pts)) ; 
H_freq_12 = (fft(h(2,:),N_pts)) ;
H_freq_21 = (fft(h(3,:),N_pts)) ; 
H_freq_22 = (fft(h(4,:),N_pts)) ; 



%%

tic;
% Loop over N symbol messages
for idx = 1:length(SNR) 
    
    for p = 1:1: 1*(1/N)*length(dataSymbolsIn)  
   

    % Calculate SNR from EbNo for each independent transmission link
    snrIndB = SNR(idx);
    snrLinear = 10^(0.1*snrIndB);


    % Create random bit vector to modulate
    msg =  dataSymbolsIn_mat(:,:,p) ; 

    OFDM_1_2 = [ msg'; msg'] ; %for simplicity ofdm streams are the same
    % Modulate data
    txSig = qammod(msg,modOrd,0);

    x = ifft(txSig,N_pts); 
    x_prefix =   x(end-(mu-1):1:end);
    x_input =  [x_prefix  ; x] ;   %dont reverse ordering
    
    y_rx_1 = conv(x_input,h(1,:).','same') + conv(x_input,h(3,:).','same') ; 
    
    y_rx_2 = conv(x_input,h(2,:).','same') + conv(x_input,h(4,:).','same') ; 
  
       
    output = transpose([ y_rx_1(mu+1:end) ,  y_rx_2(mu+1:end) ] ) ; 
    
    sig_power = 2*var(x_input,1); 
    %n = sqrt(sig_power)*10^(-snrIndB/20)*(randn(size(output)) + j*randn(size(output)))/sqrt(2) ;
    y_noisy = awgn(output,snrIndB,sig_power); 
    %y_noisy = output + n  ; 
    %
    
    Y_fft_noise_1 = fft(y_noisy(1,:),N_pts);
    Y_fft_noise_2 = fft(y_noisy(2,:),N_pts);
    
    Y_fft_noisy = [ Y_fft_noise_1 ; Y_fft_noise_2  ] ; 
       
    % Add noise to faded data   y = H* x + n
    %n = sqrt(sig_power*10^(-snrIndB/10))*(randn(N_pts, 1) + j*randn(N_pts, 1))/sqrt(2) ; 
    %y_n  = H*x_input +  n ;
    
    %output = H*x_input;
     
    Mat_64 = reshape([ H_freq_11 ;H_freq_21 ;H_freq_12 ; H_freq_22 ] , [M_r,M_t,N_pts]) ; 
    
    for kk = 1:N_pts 
          
         H = Mat_64(:,:,kk) ; 
         
         y = Y_fft_noisy(:,kk);%y_noisy(:,kk); 
         
         % no scheme
         rx_BASELINE = qamdemod(y,modOrd); 
         BER_BASELINE(idx,p, kk) = biterr(rx_BASELINE,OFDM_1_2(:,kk));
         
         %ZERO FORCING    
         W_ZF = pinv(H) ; 
         %Demodulate data
         rxSig_ZF = qamdemod(W_ZF*y,modOrd);
         BER_ZF(idx,p,kk)  = biterr(rxSig_ZF,OFDM_1_2(:,kk));
         
         
         %MMSE 
         N_0 = sig_power*10^(-snrIndB/10); 
         W_MMSE = pinv(H'*H + N_0*eye(size(H'*H )))*H' ; 
         rxSig_MMSE = qamdemod(W_MMSE*y,modOrd);
         BER_MMSE(idx,p,kk)  = biterr(rxSig_MMSE,OFDM_1_2(:,kk)) ;  
         
         
    end
    BER_BASELINE = squeeze(sum( BER_BASELINE, 3));    
    BER_ZF = squeeze(sum ( BER_ZF, 3));
    BER_MMSE = squeeze(sum ( BER_MMSE, 3));
    
    
    
     
     
    end
end    


%%

% Set up a figure for visualizing BER results

f1 = figure; %('Visible','off');
semilogy(SNR(:), mean(BER_ZF,2)/(2*N*k), 'r*' , ...
       SNR(:),   mean(BER_MMSE,2)/(2*N*k), 'b+' , ...
       SNR(:),   mean(BER_BASELINE,2)/(2*N*k), 'k^' );
xlabel('SNR (dB)');
ylabel('AVERAGE BER');
title(['OFDM WITH MIMO',' ', num2str(modOrd),'QAM',' ','N =',num2str(N_pts),...
             ' ',' ' , ' \mu = ',num2str(mu), ' ']);  
legend('ZF', 'MMSE','NO EQ','Location','best');
grid on;

saveas(f1,strcat('OFDM_N_MIMO', num2str(modOrd),'QAM_','N',num2str(N_pts),...
              '_mu',num2str(mu)),'png');


toc;

mu_ZF =   mean(BER_ZF,2)/(2*N*k) ; 
mu_MMSE =   mean(BER_MMSE,2)/(2*N*k) ; 
mu_NO_EQ =   mean(BER_BASELINE,2)/(2*N*k) ; 




%%




end

