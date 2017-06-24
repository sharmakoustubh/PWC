close all;clear all;clc;

TransTime=8;
M=78;
N=93;
DataLength=7254;

[error,OutData,crc_result] = Receiving(DataLength,TransTime);

row=78;
col=93;
pic=reshape(OutData,row,col);
imshow(pic)