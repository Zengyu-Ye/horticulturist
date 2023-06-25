%[Au, Fs]=audioread('C:\Users\CDQ\Desktop\output.mp3');   % Fs ������ 44100
Fs=100;
Au = amplitudeA(200:600,1);
[B, F, T, P] = spectrogram(Au,32,16,32,Fs);   % B��F��С��T��С�е�Ƶ�ʷ�ֵ��P�Ƕ�Ӧ���������ܶ�
figure
imagesc(T,F,C);
set(gca,'YDir','normal')
colorbar;
xlabel('ʱ�� t/s');
ylabel('Ƶ�� f/Hz');
title('��ʱ����ҶʱƵͼ');