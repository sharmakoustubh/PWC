function Tx_func(message)

fs = 44100;               
fc = 10e03;               
Nsc = 128;                
Ncp = 20;                 
N = Nsc+Ncp;              
Tsym = 58e-3;             
Rs = 1/(Tsym/Nsc);      


%% CRC
CRC_bits = comm.CRCGenerator([16 15 2 0],'CheckSumsPerFrame',1);
CRC_data = step(CRC_bits, message')';
Data_len = numel(CRC_data);

%% Convolutional coding
Trellis = poly2trellis(6,[77 45]);
Coded_bits = convenc(CRC_data,Trellis);

%% QPSK 
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

%% add Pilots 
x = zeros(1,Nsc);
randn('state',100);
P = sign(randn(1,Nsc/2));
x(1:2:end) = 2*P;

Pilot_QPSK_symbols = [x QPSK_symbols];                                                  


if mod(Data_len,Nsc) == 0
    Zero_pad = 0;
else
    Zero_pad = Nsc - mod(Data_len,Nsc);
end

Pilot_Qpsksymbols_Zeropadding = [Pilot_QPSK_symbols zeros(1,Zero_pad)];

columns = numel(Pilot_Qpsksymbols_Zeropadding)/Nsc;
Data_before_IFFT1 = reshape(Pilot_Qpsksymbols_Zeropadding,Nsc,columns);
Data_before_IFFT2 = transpose(Data_before_IFFT1); 

Data_after_IFFT  = zeros(columns,Nsc);
for i= 1:columns
    Data_after_IFFT(i,:) = ifft(Data_before_IFFT2(i,:));
end

Data_after_IFFT_stream = reshape(transpose(Data_after_IFFT),1,numel(Data_after_IFFT)); 

%% add CP
Data_withCP = [];

for i=1:Nsc:numel(Data_after_IFFT_stream)
   Data_withCP1 = Data_after_IFFT_stream(i:i+Nsc-1);
   Data_withCP2 = [Data_withCP1(Nsc-Ncp+1:end) Data_withCP1];
   Data_withCP = [Data_withCP Data_withCP2];
end

%% D/A Converter
Tx_data = interp(Data_withCP,round(fs/Rs));

%% UP-conversion
t = (0:length(Tx_data)-1)/fs;                                         
Tx_signal = sqrt(2)*real(Tx_data.*exp(1i*2*pi*fc*t));                 

%% audioplayer

y=audioplayer(Tx_signal,fs);
playblocking(y);

end
