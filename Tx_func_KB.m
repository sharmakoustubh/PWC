function Tx_func(input1)

fs = 44100;                % signal frequency
fc = 10e03;                % carrier frequency
Nsc = 128;                 % number of subcarriers
Ncp = 20;                  % length of cyclic prefix
N = Nsc+Ncp;              % symbol length              
Tsym = 58e-3;             % 1 ofdm symbol time
Rs = 1/(Tsym/Nsc);      % ofdm signal sampling rate


%% add [CRC] to picture
CRC_bits = comm.CRCGenerator([16 15 2 0],'CheckSumsPerFrame',1);
CRC_data = step(CRC_bits, input1')';
Data_len = numel(CRC_data);                                                    % used in STAGE 4 (QPSK-Mapping)

%% --- STAGE 3: add (77,45) R=1/2 Concolutional Encoder
Trellis = poly2trellis(6,[77 45]);
Coded_bits = convenc(CRC_data,Trellis);

%% --- STAGE 4: add [QPSK] Mapping 
Data_bits = reshape(Coded_bits,2,Data_len);
QPSK_symbols = zeros(1,Data_len);

for i=1:Data_len
    if(Data_bits(:,i) == [0;0])
        QPSK_symbols(1,i) = 1+1i;
    end
    if(Data_bits(:,i) == [0;1])
        QPSK_symbols(1,i) = 1-1i;
    end
    if (Data_bits(:,i) == [1;0])
        QPSK_symbols(1,i) = -1+1i;
    end
    if (Data_bits(:,i) == [1;1])
        QPSK_symbols(1,i) = -1-1i;
    end
end

%% add Pilots (Overhead) 
x = zeros(1,Nsc);
randn('state',100);
P = sign(randn(1,Nsc/2));
x(1:2:end) = 2*P;

Pilot_QPSK_symbols = [x QPSK_symbols];                                                     % [pilot qpsk-signal]

%  conversion
if mod(Data_len,Nsc) == 0
    Zero_pad = 0;
else
    Zero_pad = Nsc - mod(Data_len,Nsc);
end

Pilot_Qpsksymbols_Zeropadding = [Pilot_QPSK_symbols zeros(1,Zero_pad)];

columns = numel(Pilot_Qpsksymbols_Zeropadding)/Nsc;
Data_before_IFFT1 = reshape(Pilot_Qpsksymbols_Zeropadding,Nsc,columns);
Data_before_IFFT2 = transpose(Data_before_IFFT1); %transp(reshape(Qpsk_pilot_zeros,N_sc,N_ofdm)); kb change 34 june

% --- STAGE 5c: [IFFT] 
Data_after_IFFT  = zeros(columns,Nsc);
for i= 1:columns
    Data_after_IFFT(i,:) = ifft(Data_before_IFFT2(i,:));
end

% --- STAGE 5d: [P/S] conversion
Data_after_IFFT_stream = reshape(transpose(Data_after_IFFT),1,numel(Data_after_IFFT)); % reshape(transp(Qpsk_ifft),1,numel(Qpsk_ifft));Kb change 24 june 

%% --- STAGE 6: add [CP]
Data_withCP = [];

for i=1:Nsc:numel(Data_after_IFFT_stream)
   Data_withCP1 = Data_after_IFFT_stream(i:i+Nsc-1);
   Data_withCP2 = [Data_withCP1(Nsc-Ncp+1:end) Data_withCP1];
   Data_withCP = [Data_withCP Data_withCP2];
end

%% --- STAGE 7: [D/A] Conversion
Tx_data = interp(Data_withCP,round(fs/Rs));

%% --- STAGE 8: UP-conversion
t = (0:length(Tx_data)-1)/fs;                                         % time vector for UP-conversion
Tx_signal = sqrt(2)*real(Tx_data.*exp(1i*2*pi*fc*t));                    % final transmitted signal

%% --- STAGE 9: Wavplay
% wavplay(Tx_sig,F_s); ; % kb change 24 june
y=audioplayer(Tx_signal,fs);
playblocking(y);

end
