fs=100; %Sampling rate
wp=1.5; %ͨPass-band frequency
ws=4;   %Stop-band frequency
rp=1;   %Max fading on pass-band
as=80;  %Min fading on stop-band


%Butterworth filtering
amplitude_processed = wifi_butterworth(amplitudeA,fs,wp,ws,rp,as);
amplitude_pca = wifi_butterworth(out,fs,wp,ws,rp,as);
avg_amplitude = mean(amplitude_processed(300:900,10));
hold on
subplot(3,1,1);
plot(amplitude_pca(:,:))
subplot(3,1,2);
plot(amplitude_processed(:,:));

%����PCA��ά
[coeff,score,latent] = pca(amplitude_processed);
first_pca = score(:,1);
second_pca = score(:,2);
sequence = medfilt1(first_pca,20);
hold on
subplot(3,1,3);
plot(sequence(:,:))