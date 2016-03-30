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

%useful reference for different receivers (ZF, MMSE, precoding)
%http://www.mathworks.com/help/comm/examples/spatial-multiplexing.html


%OFDM 
% 802.11a OFDM symbol, 3 different channels, this time single path, but frequency selective.
% Implement zero-forcing and MMSE equalizer, and as before,

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


mimoChan = comm.MIMOChannel(...
  'SampleRate',                1000,...
  'PathDelays',                [0 1e-3],...
  'AveragePathGains',          [3 5],...
  'NormalizePathGains',        false,...
  'MaximumDopplerShift',       5,...
  'TransmitCorrelationMatrix', cat(3, eye(2), [1 0.1;0.1 1]),...
  'ReceiveCorrelationMatrix',  cat(3, [1 0.2;0.2 1], eye(2)),...
  'RandomStream',              'mt19937ar with seed',...
  'Seed',                      33,...
  'PathGainsOutputPort',       true);

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
