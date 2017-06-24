function DecData = Receiver(LengthData)   
% ENohnyaket, Sekyere and Gabriel
% Paramters for OFDM Receiver
    Nsc=128;
    Ncp=20;
    N=Nsc+Ncp; 
    Fc=10000;
    Tofdm=0.058;
    Fs = 44100;
    SRofdm=1/(Tofdm/Nsc);
    Recordtime =8; % 6 seconds record time
    %% Recording mesage for a given time Recordtime
    R=wavrecord(Recordtime*Fs,Fs,1);
    R=R';
    t = (1:length(R))/Fs;

    %% Down Converting
    IQ=R.*cos(2*pi*Fc*t).*sqrt(2);
    RQ=R.*(-sin(2*pi*Fc*t).*sqrt(2));
    X=IQ+1i*RQ; 

     %% LOWPASS FILTER
    [p,q]=butter(8,SRofdm/Fs,'low'); 
     r_filt=filter(p,q,X);

    %% Sampling
    samPsym=round(Fs/SRofdm);
    r_samp=r_filt(10:samPsym:end);
    
    %% Synchronization
    Tper=Nsc/2;
    miu=zeros(1,length(r_samp)-2*Tper);
    for T=1:(length(r_samp)-2*Tper)
        miu(T)=sum(r_samp(T:T+Tper-1).*conj(r_samp(T+Tper:T+2*Tper-1)))/(sqrt(sum(abs(r_samp(T:T+Tper-1).^2)))*sqrt(sum(abs(r_samp(T+Tper:T+2*Tper-1).^2))));
    end
     
    Pilot_position=find(miu==max(miu));
    %% Removing Cyclic Prefix
    Nofdm=floor(length(r_samp(Pilot_position-Ncp-1:end))/N);
   
    r_ofdm=zeros(Nofdm,Nsc);
    for ii=1:Nofdm
        r_ofdm(ii,:)=r_samp(1,(Pilot_position+N*(ii-1)):(Pilot_position+N*(ii-1)+Nsc-1)); %Nx128
    end
    %################### FFT ###########################################
    
    r_ofdm_fft=zeros(Nofdm,Nsc);
    for ii=1:Nofdm
        r_ofdm_fft(ii,:)=fft(r_ofdm(ii,:)); %Nx128
    end
    r_ofdm_fft_temp=transp(r_ofdm_fft); %128xN
    r_msg_fft=reshape(r_ofdm_fft_temp,1,numel(r_ofdm_fft_temp));

    pilot_fft=r_msg_fft(1:Nsc);
    r_msg=r_msg_fft(Nsc+1:end);

    %% Pilot  Channel estimation
    %% Pilot Generation
    randn('state',100);
    P=sign(randn(1,Nsc/2)); %1x64
    %% Channel Estimation
    Chan=pilot_fft(1:2:end)./(2*P);    %1x64./1x64=1x64
    ChanInter=interp1((1:2:Nsc),Chan,(1:Nsc)); %1-D interpolation
    ChanInter(128)=Chan(64); 

    for i=1:Nofdm-1
        r_qpsk_afterCH(Nsc*(i-1)+1:Nsc*i)=r_msg(Nsc*(i-1)+1:Nsc*i)./ChanInter; 
    end
    %% QPSK demodulation
    r_qpsk_temp=r_qpsk_afterCH;
    zI=real(r_qpsk_temp);
    zQ=imag(r_qpsk_temp);
    zI(zI>0)=0;
    zI(zI<0)=1;
    zQ(zQ>0)=0;
    zQ(zQ<0)=1;
    r_qpsk_de=[zI;zQ];
    r_VA=reshape(r_qpsk_de,1,numel([zI,zQ]));
    %% Viterbi Decoding
    trellis=poly2trellis(6,[77 45]);
    tblen=6*5;
    r_ascii=vitdec(r_VA,trellis,tblen,'term','hard');
    r_ascii_fix=r_ascii(1:LengthData+16);

    %% CRC decoding
    Hdet = comm.CRCDetector([16 15 2 0],'CheckSumsPerFrame',1);
    [TempData,CRCcheck]=step(Hdet,r_ascii_fix');
    DecData=TempData';        
end