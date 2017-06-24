function Transmitter(SendData)

%Ardiana,Jeena & Raji

%% Variables
Fs=44100; % sampling frequency
fc= 5000; % carrier frequency
Nframes= 10; %Number of frames to be transmitted
N_sc=128;  %Number of subcarriers
N_cp=20;  % Number of cyclic prefix
Tsymbol=0.058; % 1 OFDM symbol period
T_sc=Tsymbol/N_sc; % period of 1 subcarrier
deltaf=1/Tsymbol;%sampling rate or Bandwidth
Rs_ofdm=deltaf*N_sc;
%% CRC Genrator
H = comm.CRCGenerator([16 15 2 0],'CheckSumsPerFrame',1); % To append 16 bits
CRC_tx=step(H,SendData')';
SendData_lenght=numel(CRC_tx); 

%% CONVOLUTIONAL ENCODER
trellis=poly2trellis(6,[77,45]); % 6= 5 delays + 1 i/p, 77&45 means 2 o/p
coded_bits=convenc(CRC_tx,trellis);% 77=111 111 means all 5 bits affect 1st o/p, 45=100 101 means 2 delay dont affect 2 o/p

%% Zero Padding
length_zero=rem(length(coded_bits),N_sc*2);%length of zeros to make nultiple of N_sc
num_zero=2*N_sc-length_zero;%number of zeros to be added
zero_padd=zeros(1,round(num_zero));%% generation of zeros
coded_bits=[coded_bits zero_padd];% zero padding/appending

%% QPSK MODULATION
qpsk_bits=zeros(1,length(coded_bits)/2);
for k=1:length(qpsk_bits)
    if(coded_bits(2*k-1)==0&&coded_bits(2*k)==0)
        qpsk_bits(k)=1+1i;   
    end
    if(coded_bits(2*k-1)==0&&coded_bits(2*k)==1)
        qpsk_bits(k)=1-1i;  
    end
    if(coded_bits(2*k-1)==1&&coded_bits(2*k)==0)
        qpsk_bits(k)=-1+1i ;  
    end
    if(coded_bits(2*k-1)==1&&coded_bits(2*k)==1)
        qpsk_bits(k)=-1-1i;
    end         
end 
QPSK=qpsk_bits;

%% SERIAL TO PARALLEL CONVERSION
% Partitioning of message to blocks of N_sc symbols
% Each subcarrie is modulated as if it was an individual channel
size_message=length(QPSK); % 224128
block_size=floor(size_message/N_sc) %  blocks of msg/channels --> 1751
% Array of zeros with 128 rows and 1751 columns
QPSK_ifft=zeros(N_sc,block_size);
for n=1:block_size % serial to parallel and IFFT
    Ser2Par_signal=QPSK((1:N_sc)+ N_sc *(n-1)); 
    % modulate each subchannel onto the appropiate carrier
    QPSK_ifft(:,n)=ifft(Ser2Par_signal.');  
end

%% CP ADDITION
% A cyclic prefix a repetition of the first section of a symbol that is
% appended to the end of the symbol
OFDMmsg=zeros(N_sc+N_cp,block_size); % 148x1751
for hh=1:block_size % loop through the block size "individual channels"
    OFDMmsg(:,hh)=[QPSK_ifft(end-N_cp+1:end,hh); QPSK_ifft(:,hh)]; % CP+Sym
end
% figure(3)
% plot(real(OFDMmsg))
% title('OFDMmsg with CP Serial/Parallel');grid on;

%% PARALLEL TO SERIAL CONVERSION
% Suming all sub carriers and combining them into one signal
Par2ser_signal=reshape(OFDMmsg,1,[]);%% ADDITION OF PILOTS
%% GENERATION OF PILOTS
randn('state',100);
P=sign(randn(1,N_sc/2)); %1x64
x=zeros(1,N_sc);
u=2*P;
x(1:2:end)=2*P; %Every second symbol equals 0 (1x128 pilot symbols)
pilot_ifft=ifft(x); % 1x128
CP_pilot=pilot_ifft(end-N_cp+1:end);
pilot_CP=[CP_pilot pilot_ifft];
%figure(5);plot(real(pilot_CP));title('PilotSignal');grid on;
msg_plus_pilot=[ pilot_CP  Par2ser_signal]; % appending qpsk bits after the pilots
%figure(6);plot(abs(msg_plus_pilot));title('OFDMmsg plus PilotSignal');grid on;msg,1,[]);%% serial message --> 148x1751
%figure(4);plot(real(Par2ser_signal));title('OFDMmsg with CP');grid on;
%% DIGITAL TO ANALOG CONVERSION
%samples/sec*symbols/sec= samples/OFDM symbol
Nsamples= round(Fs/Rs_ofdm); % greatest integer
Interpolated_signal= interp(msg_plus_pilot,Nsamples);
[A,B]=butter(8,(Rs_ofdm/Fs));
%  DAC_signal=Interpolated_signal;
 DAC_signal=filter(A,B,Interpolated_signal);
%figure(7);plot(abs(DAC_signal));title('DAC Signal');grid on;

%% UP CONVERSION
% t=0:T_sc:T_sc*length(DAC_signal)-T_sc;
t=(0:length(DAC_signal)-1)/Fs;
I=cos(2*pi*fc*t).*sqrt(2);   %I carrier
Q=sin(2*pi*fc*t).*sqrt(2);   %Q carrier
BI=real(DAC_signal).*I;
BQ=imag(DAC_signal).*(-Q);
%OFDM_signal=BI+1i*BQ;
OFDM_signal=BI+1i*BQ;
% figure(8);plot(t,real(OFDM_signal)); xlabel('Time'); ylabel('Amplitude');title('OFDM Signal');grid on;
% figure(9);psd(real(OFDM_signal));title('POWER SPECTRUM');grid on;
R=real(OFDM_signal);

%% Transmitting
wavplay(R,Fs);
end