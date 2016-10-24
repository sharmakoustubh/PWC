clc;
close all;
clear all;
format long;

%Ardiana Osmani & Sujitha Krishnamoorthy
%% Task 1
fs= 44100; % sampling frequency 
fc=4000;    % Carrier frequency 
Ts= 0.0023;  %Symbol time 
N = round(Ts*fs);  % rounded nr of samples should the base pulse consist of wich is 101 101 times
Rs=fs/N ; % symbol rate= 436.63 symbols/sec
%% Generating base pulse
%
t1 = 0:Tsamp:Ts-Tsamp;
basepulse = sin(pi*t1/Ts);

% Energy of base pulse
E = 0;
for n=1:size(basepulse,2)
    E = E + (basepulse(n))^2;
end
E = E/Fs;
p = basepulse./sqrt(E); % Normalizing the base pulse

% figure
% subplot(2,1,1);
% stem(t1,p);
% title('Base pulse');
