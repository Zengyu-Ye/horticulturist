clear;
clc;
csi_trace=read_bf_file('./csi/isabelle/home/0126/apple4.dat');
L= length(csi_trace);
phaseA=zeros(L,30);
phaseA_a=zeros(L,30);
phaser_r=zeros(L,30);
phaser_a=zeros(L,30);
phaser_p=zeros(10,30);
phaser_pa=zeros(10,30);
phaseB=zeros(L,30);
phaseC=zeros(L,30);
phaseA2=zeros(L,30);
phaseB2=zeros(L,30);
phaseC2=zeros(L,30);
phase_var=zeros(30,1);
for m=1:L
    csi_entry=csi_trace{m};
    csi=get_scaled_csi(csi_entry);
    phaseA_r = squeeze(csi(1,1,:));
    %phaseB = squeeze(csi(2,1,:));
    phaseA2_r = squeeze(csi(1,2,:));
    phaseB2 = squeeze(csi(2,2,:));
    phaseC = squeeze(csi(1,3,:));
    phaseC2 = squeeze(csi(2,3,:));
    phaser = phaseA_r-phaseC;
    
    phaseA_r = phaseA_r.';
    phaseA2_r = phaseA2_r.';
    csi_angle=angle(squeeze(csi));    %phase
    csi_amp = abs(squeeze(csi));
    phaseA(m,:)=angle(phaseA_r);
    phaseB(m,:)=angle(phaseA2_r);
    phase_diff(m,:)=phaseA(m,:)-phaseB(m,:);
    phaseA_a(m,:)=csi_amp(1,1,:);
    

   
    
    phaser_r(m,:)=angle(phaser);
    phaser_a(m,:)=abs(phaser);
    phaser_r(1000,30)=0;
    phaser_a(1000,30)=0;
    for j=1:10
        phaser_p(j,:)=abs(mean(phaser_r(100*(j-1)+1:100*j,:)));
        phaser_pa(j,:)=mean(phaser_a(100*(j-1)+1:100*j,:));
    end
    
    for i =1:30
    phase_var(i,1) = std (phaser_p(:,i));
    end

%     phaseB(m,:)=csi2(1,2,:);
%     phaseC(m,:)=csi2(1,3,:); 
%     phaseA2(m,:)=csi_angle(2,1,:);
%     phaseB2(m,:)=csi2(2,2,:);
%     phaseC2(m,:)=csi2(2,3,:);
end
fs=100; %Sampling rate
wp=1.5; %ͨPass-band frequency
ws=4;   %Stop-band frequency
rp=1;   %Max fading on pass-band
as=80;  %Min fading on stop-band


%Butterworth filtering
phase_processed = wifi_butterworth(phase_diff,fs,wp,ws,rp,as);

hold on
% polar(phaseA,phaseA_a);%,phaseA_a(:,1),,phaser_pa(:,1)
% polar(phaser_r,phaser_a,'r');
% hold on;
% plot(real(phaseA_r),imag(phaseA_r),'r');
% plot(real(phaser),imag(phaser),'b');
%plot(phaseA);
plot(180*phaser_p/pi,'.','Markersize',25);
%legend('Original phases','Processed phases between antenna1-2');
% plot(phase_var,'r');
% xlabel('Subcarrier Index');
% ylabel('Phase Difference Variance');
% print('Subcarrier selection','-dpng');