function Transmitter(Data)
	%%Parameters
	Ncp=20;
	Nsc=128;
	N=Nsc+Ncp;
	Tp=0.058;
	Fs=44100;
	Fc=10000;
	SRofdm=1/(Tp/Nsc); 
	
	%% CRC Genrator
	H = comm.CRCGenerator([16 15 2 0],'CheckSumsPerFrame',1); % Append 16 bits
	CRCData=step(H,Data')';
	DataLen=numel(CRCData); 
	
	%% Convolutional Encoder
	%% Rate 1/2
	trellis=poly2trellis(6,[77 45]);
	CodedBits=convenc(CRCData,trellis);
	% ==> 2DataLen;real
	
	%% QPSK Modulation

    symb = [1+1j,-1+1j,1-1j,-1-1j];%./sqrt(2);
    M = log2(length(symb));
    data_bits=reshape(CodedBits,M,DataLen);
    symbol_bits = 2.^[0:M-1]*data_bits;
    QPSK_Symbols=symb(symbol_bits+1);
	
	%% Pilot Generation
	randn('state',100);
	P=sign(randn(1,Nsc/2)); 
	x=zeros(1,Nsc);
	x(1:2:end)=2*P;
	
	%% Pilot Insertion
	PilotInsertion=[x QPSK_Symbols]; 
	
	%% IFFT
	if mod(DataLen,Nsc)==0
		Nzero=0;
	else
	Nzero=Nsc-mod((DataLen),Nsc); 
	end
	% For Example DataLen=688, mod(688,128)=48, add 80 zeros to form 6 complex OFDM symbol.
	PilotInsertions=[PilotInsertion zeros(1,Nzero)]; 
	% ==> DataLen+128+Nzero;complex
	
	OFDMLeng=numel(PilotInsertions)/Nsc; 
	PilotInsertion1=reshape(PilotInsertions,Nsc,OFDMLeng);%128xN
	PilotInsertion2=transpose(PilotInsertion1);% Nx128 
	
    PilotIFFT=zeros(OFDMLeng,Nsc); % Memory for message
	
    for i=1:OFDMLeng
        PilotIFFT(i,:)=ifft(PilotInsertion2(i,:)); %6x128 9x128
    end
	
	%% Cyclic Prefix Extension
	TXsignalCP=zeros(OFDMLeng,N); % Memory for CP data
	for i=1:OFDMLeng
		TXsignalCP(i,:)=[PilotIFFT(i,Nsc-Ncp+1:end) PilotIFFT(i,1:end)];
	end
	
	%% Reshape data back
	Transdata=reshape(transpose(TXsignalCP),1,numel(TXsignalCP));
	TXsignal=interp(Transdata,round(Fs/SRofdm));% Interpolating Data
	
	%% Up Conversion
	t=(0:length(TXsignal)-1)/Fs;
	rS=sqrt(2)*real(TXsignal.*exp(1i*2*pi*Fc*t)); 

	%% Transmitting
	% wavplay(rS,Fs);
    y=audioplayer(rS,Fs);
    playblocking(y);

end