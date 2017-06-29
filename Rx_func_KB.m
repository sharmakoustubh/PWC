function[SignalwithCRC, Error] =Rx_func_KB(message)                              
                                                                           
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

%% Down Conversion
Audio_sig = R';
t1 = (1:length(Audio_sig))/fs;

I_part_down = Audio_sig.*cos(2*pi*fc*t1).*sqrt(2);
Q_part_down = Audio_sig.*(-sin(2*pi*fc*t1).*sqrt(2));

QPSK_down = I_part_down + 1j*Q_part_down;

[B,A] = butter(8,Rs/fs,'low');                                
FilteredSignal = filter(B,A,QPSK_down);                                        

%%  Sampling
Ts = round((Tsym/Nsc)*fs);
SampledSignal = FilteredSignal(1:Ts:end);                                       

%% Channel

T_per = Nsc/2;
delay = zeros(1,length(SampledSignal)-2*T_per);
 
for t1=1:(length(SampledSignal)-Nsc)
    m=SampledSignal(t1:t1+(Nsc/2)-1);
    n=SampledSignal(t1+(Nsc/2):t1+Nsc-1);
    n1=conj(n);
    delay(t1)=abs(sum(m.*n1)/(sqrt(sum((abs(m).^2)))*sqrt(sum((abs(n1).^2)))));
end

%% Removing cyclic prefix and FFT

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

Full_FFT_Signal = reshape(transpose(FFT_Signal),1,numel(FFT_Signal)); 

%% Pilot generation
randn('state',100);
P = sign(randn(1,Nsc/2));

FFT_Pilot = Full_FFT_Signal(1:Nsc);
Channel1 = FFT_Pilot(1:2:end)./(2*P);
Channel2 = interp1((1:2:Nsc),Channel1,(1:Nsc));
Channel2(128) = Channel1(63);                                           

for i=1:OFDM_Blocks-1
    QPSK_Symbols(Nsc*(i-1)+1:Nsc*i) = Full_FFT_Signal(Nsc*(i-1)+1:Nsc*i)./Channel2;
end

%% QPSK demodulation
Signal_Transmitted = QPSK_Symbols(Nsc+1:end);                                

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


%% Convolutional Decoder
const_length=6; 
Trellis = poly2trellis(const_length, [77 45]);
TB_length = 6*5;
decoded_bits = vitdec(Demod_bits_stream, Trellis, TB_length,'term','hard');
decoded_bits_fix = decoded_bits(1:message + 16); 

%% CRC detedction
CRC_Detector = comm.CRCDetector([16 15 2 0],'CheckSumsPerFrame',1);
[Signal_CRC, Error] = step(CRC_Detector, decoded_bits_fix');
SignalwithCRC = Signal_CRC';

end