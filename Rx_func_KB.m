function[Rx_ncrc, Crc_error] =Rx_func_KB(input2)                              
                                                                           
fs = 44100;                % signal frequency
fc = 10e03;                % carrier frequency
Nsc = 128;                 % number of subcarriers
Ncp = 20;                  % length of cyclic prefix
N = Nsc+Ncp;              % symbol length              
Tsym = 58e-3;             % 1 ofdm symbol time
Rs = 1/(Tsym/Nsc);      % ofdm signal sampling rate

recorder=audiorecorder(fs,16,1);
recordblocking(recorder,18);
R=getaudiodata(recorder); 

%figure(1)
%plot(R)

%% DOWN-conversion

Rsig = R';
t1 = (1:length(Rsig))/fs;

I_part_down = Rsig.*cos(2*pi*fc*t1).*sqrt(2);
Q_part_down = Rsig.*(-sin(2*pi*fc*t1).*sqrt(2));

Rx_sig = Rx_sig_i + 1j*Rx_sig_q;

[B,A] = butter(8,Rs/fs,'low');                                
FilteredSignal = filter(B,A,Rx_sig);                                        

%%  [A/D] conversion - sampling
%N_samp = round(fs/Rs);
Ts = round((Tsym/Nsc)*fs)
SampledSignal = FilteredSignal(1:Ts:end);                                       

%% STAGE 7a': Channel estimation
T_per = Nsc/2;
delay = zeros(1,length(SampledSignal)-2*T_per);
 
for t1=1:(length(SampledSignal)-Nsc)
    
     a = SampledSignal(t1:t1+(Nsc/2));
     b = conj(SampledSignal(t1+(Nsc/2):t1+2*(Nsc/2)-1));
     c = sqrt(sum(abs(SampledSignal(t1:t1+(Nsc/2)-1).^2)));
     d = sqrt(sum(abs(SampledSignal(t1+(Nsc/2):t1+2*(Nsc/2)-1).^2)));

     nom = sum(a.*b);
     denom = c*d;
     delay(t1) = nom/denom;
end

Pilot_pos = find(delay == max(delay));
 
N_ofdm = floor(length(SampledSignal(Pilot_pos:end))/N);

%%  Removing [CP], S/P
Rx_qpsk_cp = zeros(N_ofdm,Nsc);

for i=1:N_ofdm-1
    Rx_qpsk_cp(i,:) = SampledSignal(1,(Pilot_pos + N*(i-1)):(Pilot_pos + N*(i-1) + Nsc - 1));
end

%% --- STAGE 5c': FFT
Rx_ofdm_fft = zeros(N_ofdm,Nsc);

for i=1:N_ofdm
    Rx_ofdm_fft(i,:) = fft(Rx_qpsk_cp(i,:));
end

% --- STAGE 5b': P/S
Rx_fft_ser = reshape(transpose(Rx_ofdm_fft),1,numel(Rx_ofdm_fft)); % reshape(transp(Rx_ofdm_fft),1,numel(Rx_ofdm_fft)); kb 24 june

% --- STAGE 5a': Channel Estimation based on pilots
randn('state',100);
P = sign(randn(1,Nsc/2));

Qpsk_ch_tmp = Rx_fft_ser(1:Nsc);
Alpha = Qpsk_ch_tmp(1:2:end)./(2*P);
Alpha_int = interp1((1:2:Nsc),Alpha,(1:Nsc));
Alpha_int(128) = Alpha(63);                                                % interpolation error - shift picture ->(*)

for i=1:N_ofdm-1
    Qpsk_Rx_ch(Nsc*(i-1)+1:Nsc*i) = Rx_fft_ser(Nsc*(i-1)+1:Nsc*i)./Alpha_int;
end

%% --- Stage 4': QPSK demodulation
Qpsk_Rx_tmp = Qpsk_Rx_ch(Nsc+1:end);                                      % removing interpolation shift error (*)
z_i = real(Qpsk_Rx_tmp);
z_q = imag(Qpsk_Rx_tmp);

z_i(z_i > 0) = 0;
z_i(z_i < 0) = 1;

z_q(z_q > 0) = 0;
z_q(z_q < 0) = 1;

Qpsk_Rx_demod = [z_i;z_q];
Rx_demod = reshape(Qpsk_Rx_demod,1,numel([z_i,z_q]));

%% --- STAGE 3': (77,45) R=1/2 Concolutional Decoder
code = poly2trellis(6, [77 45]);
L_tb = 6*5;
Pic_Rx_ascii = vitdec(Rx_demod, code, L_tb,'term','hard');
Pic_Rx_ascii = Pic_Rx_ascii(1:input2 + 16); 

%% --- STAGE 2': CRC detedction
Crc_detec = comm.CRCDetector([16 15 2 0],'CheckSumsPerFrame',1);
[Rx_crc, Crc_error] = step(Crc_detec, Pic_Rx_ascii');
Rx_ncrc = Rx_crc';

%% --- STAGE 1': reshape msg, show picture
% Picture_Rx_crc = reshape(Rx_ncrc,143,144);
% imshow(Picture_Rx_crc)
% disp(bi2de(Rx_ncrc(1:10)))
% disp(bi2de(Rx_ncrc(11:end)))
% disp('Number of CRC errors:')
% disp(Crc_error)

end % function ending here