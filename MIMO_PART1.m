function [ BERR_ZF,BERR_MMSE,BERR_PRECODING,BERR_BASELINE ,BER_BASELINE] = MIMO_PART1( H, SNR,str_type,modOrd)

%http://www.mathworks.com/help/comm/examples/spatial-multiplexing.html
% channel side info 
%generate QPSK data
%M                  % Number of transmit antennas
  %N                % Number of receive antennas
[M N ] = size(H);
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
n = k*4*1000;                  % Number of bits to process
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
%SNR = EbNoVec + 10*log10(modOrd);

tic;
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

          
    % Add noise to faded data   y = H* x + n
    n = 10^(-snrIndB/20)*(randn(M, 1) + j*randn(M, 1))/sqrt(2) ; 
    y  = H*txSig +  n ;
     
    rx_BASELINE = qamdemod(y,modOrd); 
    BER_BASELINE(idx,p) = biterr(rx_BASELINE,msg);
      
    %PRECODING (SVD)
    [U, SIGMA, V] = svd(H) ;
    x_tilde = V* txSig;
       
       
    y_tilde = SIGMA* x_tilde +  U'*n;
       
    %Demodulate data
    rxSig_PRECODING = qamdemod(U*y_tilde,modOrd);
    %BER
    BER_PRECODING(idx,p)  = biterr(rxSig_PRECODING,msg);
   
    
       
     
    %ZERO FORCING    
    W_ZF = pinv(H) ; 
    %Demodulate data
    rxSig_ZF = qamdemod(W_ZF*y,modOrd);
    BER_ZF(idx,p)  = biterr(rxSig_ZF,msg);
    
       
    %MMSE 
    N_0 = 1; 
    W_MMSE = inv(H'*H + N_0)*H' ; 
    rxSig_MMSE = qamdemod(W_MMSE*y,modOrd);
    BER_MMSE(idx,p)  = biterr(rxSig_MMSE,msg) ;  
       
     
    
    end    
end    

toc;
    %%
       
% Set up a figure for visualizing BER results
f1 = figure('Visible','off');
semilogy(SNR(:), mean(BER_ZF,2)/modOrd, 'r*' , ...
       SNR(:), mean(BER_MMSE,2)/modOrd, 'b+' , ...
       SNR(:), mean(BER_PRECODING,2)/modOrd, 'mo', ...
       SNR(:), mean(BER_BASELINE,2)/modOrd, 'k^');

xlim([SNR(1)-0.01, SNR(end)]);
xlabel('SNR (dB)');
ylabel('AVERAGE BER');
title(['2x2','', ' ' , num2str(modOrd),'QAM',' System',' ' , str_type]);  
legend('ZF', 'MMSE', 'PRECODING','NO SCHEME','Location','best');
grid on;

saveas(f1,strcat('2x2','_MIMO_', num2str(modOrd),'QAM', str_type),'png');
%%



BERR_ZF = mean(BER_ZF,2)/modOrd ;
BERR_MMSE = mean(BER_MMSE,2)/modOrd;
BERR_PRECODING = mean(BER_PRECODING,2)/modOrd;
BERR_BASELINE = mean(BER_BASELINE,2)/modOrd;
      


end

