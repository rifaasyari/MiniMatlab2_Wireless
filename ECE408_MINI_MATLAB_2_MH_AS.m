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
 

%http://www.mathworks.com/help/comm/examples/spatial-multiplexing.html
% channel side info 
%generate QPSK data
N = 2;                  % Number of transmit antennas
M = 2;                  % Number of receive antennas
EbNoVec = 2:3:8;        % Eb/No in dB
modOrd = 2;             % constellation size = 2^modOrd

a local random stream to be used by random number generators for
% repeatability.
hStr = RandStream('mt19937ar');

% Create PSK modulator and demodulator System objects
hMod   = comm.PSKModulator(...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitInput',         true);
hDemod = comm.PSKDemodulator( ...
            'ModulationOrder',  2^modOrd, ...
            'PhaseOffset',      0, ...
            'BitOutput',        true);

% Create error rate calculation System objects for 3 different receivers
hZFBERCalc   = comm.ErrorRate;
hMMSEBERCalc = comm.ErrorRate;
hMLBERCalc   = comm.ErrorRate;

% Get all bit and symbol combinations for ML receiver
allBits  = de2bi(0:2^(modOrd*N)-1, 'left-msb')';
allTxSig = reshape(step(hMod, allBits(:)), N, 2^(modOrd*N));

% Pre-allocate variables to store BER results for speed
[BER_ZF, BER_MMSE, BER_ML] = deal(zeros(length(EbNoVec), 3));



% Set up a figure for visualizing BER results
h = gcf;
grid on;
hold on;
ax = gca;
ax.YScale = 'log'
xlim([EbNoVec(1)-0.01 EbNoVec(end)]);
ylim([1e-3 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
h.NumberTitle = 'off';
h.Renderer = 'zbuffer';
h.Name = 'Spatial Multiplexing';
title('2x2 Uncoded QPSK System');

% Loop over selected EbNo points
for idx = 1:length(EbNoVec)
    % Reset error rate calculation System objects
    reset(hZFBERCalc);
    reset(hMMSEBERCalc);
    reset(hMLBERCalc);

    % Calculate SNR from EbNo for each independent transmission link
    snrIndB = EbNoVec(idx) + 10*log10(modOrd);
    snrLinear = 10^(0.1*snrIndB);

    while (BER_ZF(idx, 3) < 1e5) && ((BER_MMSE(idx, 2) < 100) || ...
          (BER_ZF(idx, 2) < 100) ||  (BER_ML(idx, 2)   < 100))
        % Create random bit vector to modulate
        msg = randi(hStr, [0 1], [N*modOrd, 1]);

        % Modulate data
        txSig = step(hMod, msg);

        % Flat Rayleigh fading channel with independent links
        rayleighChan = (randn(hStr, M, N) +  1i*randn(hStr, M, N))/sqrt(2);

        % Add noise to faded data
        rxSig = awgn(rayleighChan*txSig, snrIndB, 0, hStr);

        % ZF-SIC receiver
        r = rxSig;
        H = rayleighChan; % Assume perfect channel estimation
        % Initialization
        estZF = zeros(N*modOrd, 1);
        orderVec = 1:N;
        k = N+1;
        % Start ZF nulling loop
        for n = 1:N
            % Shrink H to remove the effect of the last decoded symbol
            H = H(:, [1:k-1,k+1:end]);
            % Shrink order vector correspondingly
            orderVec = orderVec(1, [1:k-1,k+1:end]);
            % Select the next symbol to be decoded
            G = (H'*H) \ eye(N-n+1); % Same as inv(H'*H), but faster
            [~, k] = min(diag(G));
            symNum = orderVec(k);

            % Hard decode the selected symbol
            decBits = step(hDemod, G(k,:) * H' * r);
            estZF(modOrd * (symNum-1) + (1:modOrd)) = decBits;

            % Subtract the effect of the last decoded symbol from r
            if n < N
                r = r - H(:, k) * step(hMod, decBits);
            end
        end

        % MMSE-SIC receiver
        r = rxSig;
        H = rayleighChan;
        % Initialization
        estMMSE = zeros(N*modOrd, 1);
        orderVec = 1:N;
        k = N+1;
        % Start MMSE nulling loop
        for n = 1:N
            H = H(:, [1:k-1,k+1:end]);
            orderVec = orderVec(1, [1:k-1,k+1:end]);
            % Order algorithm (matrix G calculation) is the only difference
            % with the ZF-SIC receiver
            G = (H'*H + ((N-n+1)/snrLinear)*eye(N-n+1)) \ eye(N-n+1);
            [~, k] = min(diag(G));
            symNum = orderVec(k);

            decBits = step(hDemod, G(k,:) * H' * r);
            estMMSE(modOrd * (symNum-1) + (1:modOrd)) = decBits;

            if n < N
                r = r - H(:, k) * step(hMod, decBits);
            end
        end

        % ML receiver
        r = rxSig;
        H = rayleighChan;
        [~, k] = min(sum(abs(repmat(r,[1,2^(modOrd*N)]) - H*allTxSig).^2));
        estML = allBits(:,k);

        % Update BER
        BER_ZF(  idx, :) = step(hZFBERCalc,   msg, estZF);
        BER_MMSE(idx, :) = step(hMMSEBERCalc, msg, estMMSE);
        BER_ML(  idx, :) = step(hMLBERCalc,   msg, estML);
    end

    % Plot results
    semilogy(EbNoVec(1:idx), BER_ZF(  1:idx, 1), 'r*', ...
             EbNoVec(1:idx), BER_MMSE(1:idx, 1), 'bo', ...
             EbNoVec(1:idx), BER_ML(  1:idx, 1), 'gs');
    legend('ZF-SIC', 'MMSE-SIC', 'ML');
    drawnow;
end

% Draw the lines
semilogy(EbNoVec, BER_ZF(  :, 1), 'r-', ...
         EbNoVec, BER_MMSE(:, 1), 'b-', ...
         EbNoVec, BER_ML(  :, 1), 'g-');
hold off;

%still need to use precoding!!!!!


%% PART 2 - OFDM 












%% PART 3 - OFDM & MIMO

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
