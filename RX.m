close all;
clc;

% Ardiana, Jeena, Raji

%% Receive

Transmission_t=8;
Lenght=7254;

[err,FinalData,crc]=Receive(Lenght,Transmission_t);

row=78;
col=93;
pic=reshape(FinalData,row,col);
imshow(pic)