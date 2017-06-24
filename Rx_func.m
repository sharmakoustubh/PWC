function[Rx_ncrc, Crc_error] =Rx_func(input2)                              % input 2 = number of transmitted payload bits
                                                                           % without 16 CRC bits
%% --- Receiver -----------------------------------------------------------

%% ---Variable Definitions-------------------------------------------------
F_s = 44100;                % signal frequency
F_c = 10e03;                % carrier frequency
N_sc = 128;                 % number of subcarriers
L_cp = 20;                  % length of cyclic prefix
N = N_sc+L_cp;              % symbol length              
T_ofdm = 58e-3;             % 1 ofdm symbol time
R_s = 1/(T_ofdm/N_sc);      % ofdm signal sampling rate

Picture_Rx = wavrecord(F_s,F_s,1);
%figure(1)
%plot(Picture_Rx)

%% --- STAGE 8b': DOWN-conversion
Rec_sig = Picture_Rx';
t = (1:length(Rec_sig))/F_s;

C_i = cos(2*pi*F_c*t).*sqrt(2);
C_q = -sin(2*pi*F_c*t).*sqrt(2);

Rx_sig_i = Rec_sig.*C_i;
Rx_sig_q = Rec_sig.*C_q;

Rx_sig = Rx_sig_i + 1j*Rx_sig_q;

% --- STAGE 8a': Lowpass Filter                                            
[B,A] = butter(8,R_s/F_s,'low');                                
Rx_filt = filter(B,A,Rx_sig);                                              % Removing mirror freq.

%% --- STAGE 7b': [A/D] conversion - sampling
N_samp = round(F_s/R_s);
Rx_qpsk_dig = Rx_filt(1:N_samp:end);                                       % Start sampling from correct probe!!!

%% STAGE 7a': Channel estimation
T_per = N_sc/2;
miu = zeros(1,length(Rx_qpsk_dig)-2*T_per);
 
for t=1:(length(Rx_qpsk_dig)-2*T_per)
    
     a = Rx_qpsk_dig(t:t+T_per-1);
     b = conj(Rx_qpsk_dig(t+T_per:t+2*T_per-1));
     c = sqrt(sum(abs(Rx_qpsk_dig(t:t+T_per-1).^2)));
     d = sqrt(sum(abs(Rx_qpsk_dig(t+T_per:t+2*T_per-1).^2)));

     nom = sum(a.*b);
     denom = c*d;
     miu(t) = nom/denom;
end
 
Pilot_pos = find(miu == max(miu));
 
N_ofdm = floor(length(Rx_qpsk_dig(Pilot_pos:end))/N);

%% --- STAGE 6': Removing [CP], S/P
Rx_qpsk_cp = zeros(N_ofdm,N_sc);

for i=1:N_ofdm-1
    Rx_qpsk_cp(i,:) = Rx_qpsk_dig(1,(Pilot_pos + N*(i-1)):(Pilot_pos + N*(i-1) + N_sc - 1));
end

%% --- STAGE 5c': FFT
Rx_ofdm_fft = zeros(N_ofdm,N_sc);

for i=1:N_ofdm
    Rx_ofdm_fft(i,:) = fft(Rx_qpsk_cp(i,:));
end

% --- STAGE 5b': P/S
Rx_fft_ser = reshape(transp(Rx_ofdm_fft),1,numel(Rx_ofdm_fft));

% --- STAGE 5a': Channel Estimation based on pilots
randn('state',100);
P = sign(randn(1,N_sc/2));

Qpsk_ch_tmp = Rx_fft_ser(1:N_sc);
Alpha = Qpsk_ch_tmp(1:2:end)./(2*P);
Alpha_int = interp1((1:2:N_sc),Alpha,(1:N_sc));
Alpha_int(128) = Alpha(63);                                                % interpolation error - shift picture ->(*)

for i=1:N_ofdm-1
    Qpsk_Rx_ch(N_sc*(i-1)+1:N_sc*i) = Rx_fft_ser(N_sc*(i-1)+1:N_sc*i)./Alpha_int;
end

%% --- Stage 4': QPSK demodulation
Qpsk_Rx_tmp = Qpsk_Rx_ch(N_sc+1:end);                                      % removing interpolation shift error (*)
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