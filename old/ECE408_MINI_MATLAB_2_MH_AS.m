%{ 
Alexander Serrano & Max Howald
ECE 408 - WIRELESS COMMS
Prof. Keene
MiniMatlab Assignment #2 
%}

%Source: (MATHWORKS) "OFDM with MIMO Simulation"
%http://www.mathworks.com/help/comm/ug/ofdm-with-mimo-simulation.html


%MIMO 
% 2X2 , FLAT FADING GAINS, ( Pre-coding - CSIT , Zero-Forcing and MMSE - CSIR) .




%OFDM 
% 802.11a OFDM symbol, 3 different channels, this time single path, but frequency selective.
% Implement zero-forcing and MMSE equalizer, and as before,


%% PART 1 - MIMO 
clc; clear all; close force all;

warning('off','all')
% CHANNEL 1
str_type1 = 'Correlated Case, H (FULL RANK)';
SNR = 0:1:30;
H_1 =  0.1*[ 1 , 2 ; 7 , 3 ]  ;
modOrd = 4;

MIMO_PART1( H_1,SNR,str_type1,modOrd);     % 2qam

%%
MIMO_PART1( H_1,SNR,str_type1,modOrd^2);    % 4qam
MIMO_PART1( H_1,SNR,str_type1,modOrd^3);     %8 qam
MIMO_PART1( H_1,SNR,str_type1,modOrd^4);      %16 qam

%% CHANNEL 2
str_type2 = 'Uncorrelated Case, H (FULL RANK)';
H_2 = 0.1*[ 2 0 ; 0 ,3.8];
MIMO_PART1( H_2,SNR,str_type2,modOrd);     % 2qam
MIMO_PART1( H_2,SNR,str_type2,modOrd^2);    % 4qam
MIMO_PART1( H_2,SNR,str_type2,modOrd^3);     %8 qam
MIMO_PART1( H_2,SNR,str_type2,modOrd^4);      %16 qam 

%% CHANNEL 3
str_type3 = 'Correlated Case, H (NOT FULL RANK)';
H_3 = 0.1*[ 6  3; 4 ,2];
MIMO_PART1( H_3,SNR,str_type3,modOrd);     % 2qam
MIMO_PART1( H_3,SNR,str_type3,modOrd^2);    % 4qam
MIMO_PART1( H_3,SNR,str_type3,modOrd^3);     %8 qam
MIMO_PART1( H_3,SNR,str_type3,modOrd^4);      %16 qam 
    




%% PART 2 - OFDM 

clc; clear all; close force all; 
SNR = 0:1:50;
N_pts = 64;

modOrd = 4; 

mu = 16  ; % size of cyclic prefix
%h1 = [ 1 ,zeros(1,mu+1) ] ;
h2 = [ 1 , 0.1, 0.9 , zeros(1,14) ];
%h2 =  [ 1 , 0.9 , 0.1, 0, zeros(1,13) ];

%h2 = [1 ] ;
h3 = [ 0.94 , 0.17 , 0.05 0.01 0.01*ones(1,13) ];
%h1 = 0.1*[ 5, 3, 2 ,zeros(1,14) ];  %MU = length(4)

N = 64;

tic;
%2QAM
OFDM_PART2(h2, SNR,modOrd,mu,N_pts,N,'h2');

%%
OFDM_PART2(h2, SNR,modOrd,mu,N_pts,N,'h2');
OFDM_PART2(h3, SNR,modOrd,mu,N_pts,N,'h3');

toc;


%%
tic;
%4QAM
OFDM_PART2(h1, SNR,modOrd^2,mu,N_pts,N,'h1');
OFDM_PART2(h2, SNR,modOrd^2,mu,N_pts,N,'h2');
OFDM_PART2(h3, SNR,modOrd^2,mu,N_pts,N,'h3');
toc;

tic;
%8QAM
OFDM_PART2(h1, SNR,modOrd^3,mu,N_pts,N,'h1');
OFDM_PART2(h2, SNR,modOrd^3,mu,N_pts,N,'h2');
OFDM_PART2(h3, SNR,modOrd^3,mu,N_pts,N,'h3');
toc; 

tic ; 
%16QAM
OFDM_PART2(h1, SNR,modOrd^4,mu,N_pts,N,'h1');
OFDM_PART2(h2, SNR,modOrd^4,mu,N_pts,N,'h2');
OFDM_PART2(h3, SNR,modOrd^4,mu,N_pts,N,'h3');
toc;





%% PART 3 - OFDM & MIMO

%clc; clear all; close force all; 
SNR = 1:5:50  %:1:50;
N_pts = 64;

modOrd = 4; 
%
mu = 16  ; % size of cyclic prefix

M_r = 2;
M_t = 2;
h =  [ 1  0  0  0 ]  ; 
%h = [ h ; 1.05*h ; 0.975*h ; 0.99* h ]
h =  [h ;  zeros(3,4) ] ; 
N = 64;




tic;
[ mu_ZF , mu_MMSE , mu_NO_EQ ] = OFDM_N_MIMO_PART3(h, SNR,modOrd,mu,M_r,M_t,N_pts,N)
toc;

%%
OFDM_N_MIMO_PART3(h, SNR,modOrd,mu,N_pts,N)

%% OFDM modulator and demodulator in a simple, 2x2 MIMO error rate simulation

%Create an OFDM modulator and demodulator pair
qpskMod = comm.QPSKModulator;
qpskDemod = comm.QPSKDemodulator;

ofdmMod = comm.OFDMModulator('FFTLength',128,'PilotInputPort',true,...
    'PilotCarrierIndices',cat(3,[12; 40; 54; 76; 90; 118],...
    [13; 39; 55; 75; 91; 117]),'InsertDCNull',true,...
    'NumTransmitAntennas',2);
ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDemod.NumReceiveAntennas = 2;



%make AWGN channel 

ch = comm.AWGNChannel(...
    'NoiseMethod','Signal to noise ratio (Es/No)', ...
    'EsNo',30);

%dimensions of OFDM Modulator 
numData = ofdmModDim.DataInputSize(1);   % Number of data subcarriers
numSym = ofdmModDim.DataInputSize(2);    % Number of OFDM symbols
numTxAnt = ofdmModDim.DataInputSize(3);  % Number of transmit antennas


%generate data symbols
nframes = 100;
data = randi([0 3],nframes*numData,numSym,numTxAnt);

%apply modulation:
modData = step(qpskMod,data(:));
modData = reshape(modData,nframes*numData,numSym,numTxAnt);

%introduce error rate counter
err = comm.ErrorRate;



%% scenario - OFDM system over 100 frames assuming a flat, 2x2, Rayleigh fading channel.
for k = 1:nframes

    % Find row indices for kth OFDM frame
    indData = (k-1)*ofdmModDim.DataInputSize(1)+1:k*numData;

    % Generate random OFDM pilot symbols
    pilotData = complex(rand(ofdmModDim.PilotInputSize), ...
        rand(ofdmModDim.PilotInputSize));

    % Modulate QPSK symbols using OFDM
    dataOFDM = step(ofdmMod,modData(indData,:,:),pilotData);

    % Create flat, i.i.d., Rayleigh fading channel
    chGain = complex(randn(2,2),randn(2,2))/sqrt(2); % Random 2x2 channel

    % Pass OFDM signal through Rayleigh and AWGN channels
    receivedSignal = step(ch,dataOFDM*chGain);

    % Apply least squares solution to remove effects of fading channel
    rxSigMF = chGain.' \ receivedSignal.';

    % Demodulate OFDM data
    receivedOFDMData = step(ofdmDemod,rxSigMF.');

    % Demodulate QPSK data
    receivedData = step(qpskDemod,receivedOFDMData(:));

    % Compute error statistics
    dataTmp = data(indData,:,:);
    errors = step(err,dataTmp(:),receivedData);
end

fprintf('\nSymbol error rate = %d from %d errors in %d symbols\n',errors)
