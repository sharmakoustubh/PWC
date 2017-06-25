function[SignalwithCRC, Error] =Rx_func_KB(input2)                              
                                                                           
fs = 44100;                
fc = 10e03;                
Nsc = 128;                 
Ncp = 20;                  
N = Nsc+Ncp;              
Tsym = 58e-3;             
Rs = 1/(Tsym/Nsc);      

recorder=audiorecorder(fs,16,1);
recordblocking(recorder,18);
R=getaudiodata(recorder); 

%% DOWN-conversion

Rsig = R';
t1 = (1:length(Rsig))/fs;

I_part_down = Rsig.*cos(2*pi*fc*t1).*sqrt(2);
Q_part_down = Rsig.*(-sin(2*pi*fc*t1).*sqrt(2));

Rx_sig = I_part_down + 1j*Q_part_down;

[B,A] = butter(8,Rs/fs,'low');                                
FilteredSignal = filter(B,A,Rx_sig);                                        

%%  [A/D] conversion - sampling
%N_samp = round(fs/Rs);
Ts = round((Tsym/Nsc)*fs)
SampledSignal = FilteredSignal(1:Ts:end);                                       

%% STAGE 7a': Channel estimation
T_per = Nsc/2;
delay = zeros(1,length(SampledSignal)-2*T_per);
 
% for t1=1:(length(SampledSignal)-Nsc)
%     
%      m = SampledSignal(t1:t1+(Nsc/2)-1);
%      n = conj(SampledSignal(t1+(Nsc/2):t1+Nsc-1));
%      c = sqrt(sum(abs(SampledSignal(t1:t1+(Nsc/2)-1).^2)));
%      d = sqrt(sum(abs(SampledSignal(t1+(Nsc/2):t1+2*(Nsc/2)-1).^2)));
% 
%      nom = sum(m.*n);
%      denom = c*d;
%      delay(t1) = nom/denom;
% end
 
 
for t1=1:(length(SampledSignal)-Nsc)
    m=SampledSignal(t1:t1+(Nsc/2)-1);
    n=SampledSignal(t1+(Nsc/2):t1+Nsc-1);
    n1=conj(n);
    delay(t1)=abs(sum(m.*n1)/(sqrt(sum((abs(m).^2)))*sqrt(sum((abs(n1).^2)))));
end
  

% for t1=1:length(SampledSignal)-(Nsc); 
%     m=SampledSignal(t1:t1+(Nsc/2));
%     n=SampledSignal(t1+(Nsc/2):t1+Nsc);
%     n1=conj(n);
%     delay(t1)=abs(sum(m.*n1)/(sqrt(sum((abs(m).^2)))*sqrt(sum((abs(n1).^2)))));
% end

%% Removing cyclic prefix and FFT

%start_pos = find(delay == max(delay));
[Corr start_pos]= max(delay);
OFDM_Blocks = floor(length(SampledSignal(start_pos:end))/N);

OFDM_chunk = zeros(OFDM_Blocks,Nsc);

for i=1:OFDM_Blocks-1
    OFDM_chunk(i,:) = SampledSignal(1,(start_pos + N*(i-1)):(start_pos + N*(i-1) + Nsc - 1));
end

FFT_Signal = zeros(OFDM_Blocks,Nsc);

for i=1:OFDM_Blocks
    FFT_Signal(i,:) = fft(OFDM_chunk(i,:));
end

Full_FFT_Signal = reshape(transpose(FFT_Signal),1,numel(FFT_Signal)); % reshape(transp(Rx_ofdm_fft),1,numel(Rx_ofdm_fft)); kb 24 june

%% Pilot generation

randn('state',100);
P = sign(randn(1,Nsc/2));

FFT_Pilot = Full_FFT_Signal(1:Nsc);
Channel1 = FFT_Pilot(1:2:end)./(2*P);
Channel2 = interp1((1:2:Nsc),Channel1,(1:Nsc));
Channel2(128) = Channel1(63);                                                % interpolation error - shift picture ->(*)

for i=1:OFDM_Blocks-1
    QPSK_Symbols(Nsc*(i-1)+1:Nsc*i) = Full_FFT_Signal(Nsc*(i-1)+1:Nsc*i)./Channel2;
end

%% --- Stage 4': QPSK demodulation
Signal_Transmitted = QPSK_Symbols(Nsc+1:end);                                      % removing interpolation shift error (*)

I_part=sign(real(Signal_Transmitted)); 
Q_part=sign(imag(Signal_Transmitted));

for k=1:length(I_part)   
    if (I_part(k)>0)
        m(k)= 0;
    else
        m(k)=1;
    end
end
for k=1:length(Q_part)  
    if (Q_part(k)>0)
        n(k)= 0;
    else
        n(k)=1;
    end
end
    
Demod_bits = [m; n];    
Demod_bits_stream = reshape(Demod_bits,1,numel(Demod_bits));


%% (77,45) R=1/2 Concolutional Decoder
const_length=6; 
Trellis = poly2trellis(const_length, [77 45]);
TB_length = 6*5;
decoded_bits = vitdec(Demod_bits_stream, Trellis, TB_length,'term','hard');
decoded_bits_fix = decoded_bits(1:input2 + 16); 

%% CRC detedction
CRC_Detector = comm.CRCDetector([16 15 2 0],'CheckSumsPerFrame',1);
[Signal_CRC, Error] = step(CRC_Detector, decoded_bits_fix');
SignalwithCRC = Signal_CRC';

%% --- STAGE 1': reshape msg, show picture
% Picture_Rx_crc = reshape(Rx_ncrc,143,144);
% imshow(Picture_Rx_crc)
% disp(bi2de(Rx_ncrc(1:10)))
% disp(bi2de(Rx_ncrc(11:end)))
% disp('Number of CRC errors:')
% disp(Crc_error)

end % function ending here