function Tx_func(input1)
%% Wir Comm proj, final part

% ----TRANSMITER-----------------------------------------------------------
% clear all
% close all
%% ---Variable Definitions-------------------------------------------------
F_s = 44100;                % signal frequency
F_c = 10e03;                % carrier frequency
N_sc = 128;                 % number of subcarriers
L_cp = 20;                  % length of cyclic prefix
N = N_sc+L_cp;              % symbol length              
T_ofdm = 58e-3;             % 1 ofdm symbol time
R_s = 1/(T_ofdm/N_sc);      % ofdm signal sampling rate

%% --- STAGE 1: picture reading, reshaping
Pic = input1;

%% --- STAGE 2: add [CRC] to picture
Crc_gen = comm.CRCGenerator([16 15 2 0],'CheckSumsPerFrame',1);
Pic_crc = step(Crc_gen, Pic')';
L_pic = numel(Pic_crc);                                                    % used in STAGE 4 (QPSK-Mapping)

%% --- STAGE 3: add (77,45) R=1/2 Concolutional Encoder
code = poly2trellis(6,[77 45]);
Pic_crc_conv = convenc(Pic_crc,code);

%% --- STAGE 4: add [QPSK] Mapping 
Qpsk_tmp = reshape(Pic_crc_conv,2,L_pic);
Qpsk = zeros(1,L_pic);

for i=1:L_pic
    if(Qpsk_tmp(:,i) == [0;0])
        Qpsk(1,i) = 1+1i;
    end
    if(Qpsk_tmp(:,i) == [0;1])
        Qpsk(1,i) = 1-1i;
    end
    if (Qpsk_tmp(:,i) == [1;0])
        Qpsk(1,i) = -1+1i;
    end
    if (Qpsk_tmp(:,i) == [1;1])
        Qpsk(1,i) = -1-1i;
    end
end

%% --- STAGE 5a: add Pilots (Overhead) 
x = zeros(1,N_sc);
randn('state',100);
P = sign(randn(1,N_sc/2));
x(1:2:end) = 2*P;

Qpsk_pilot = [x Qpsk];                                                     % [pilot qpsk-signal]

% --- STAGE 5b: [S/P] conversion
if mod(L_pic,N_sc) == 0
    Nzero = 0;
else
    Nzero = N_sc - mod(L_pic,N_sc);
end

Qpsk_pilot_zeros = [Qpsk_pilot zeros(1,Nzero)];

N_ofdm = numel(Qpsk_pilot_zeros)/N_sc;
Qpsk_par = transpose(reshape(Qpsk_pilot_zeros,N_sc,N_ofdm)); %transp(reshape(Qpsk_pilot_zeros,N_sc,N_ofdm)); kb change 34 june

% --- STAGE 5c: [IFFT] 
Qpsk_ifft = zeros(N_ofdm,N_sc);
for i= 1:N_ofdm
    Qpsk_ifft(i,:) = ifft(Qpsk_par(i,:));
end

% --- STAGE 5d: [P/S] conversion
Qpsk_ifft_ser = reshape(transpose(Qpsk_ifft),1,numel(Qpsk_ifft)); % reshape(transp(Qpsk_ifft),1,numel(Qpsk_ifft));Kb change 24 june 

%% --- STAGE 6: add [CP]
Qpsk_cp = [];

for i=1:N_sc:numel(Qpsk_ifft_ser)
   Qpsk_cp_tmp_A = Qpsk_ifft_ser(i:i+N_sc-1);
   Qpsk_cp_tmp_B = [Qpsk_cp_tmp_A(N_sc-L_cp+1:end) Qpsk_cp_tmp_A];
   Qpsk_cp = [Qpsk_cp Qpsk_cp_tmp_B];
end

%% --- STAGE 7: [D/A] Conversion
Qpsk_cp_d2a = interp(Qpsk_cp,round(F_s/R_s));

%% --- STAGE 8: UP-conversion
t = (0:length(Qpsk_cp_d2a)-1)/F_s;                                         % time vector for UP-conversion

Tx_sig = sqrt(2)*real(Qpsk_cp_d2a.*exp(1i*2*pi*F_c*t));                    % final transmitted signal

%% --- STAGE 9: Wavplay
% wavplay(Tx_sig,F_s); ; % kb change 24 june
audioplayer(Tx_sig,F_s)
end
