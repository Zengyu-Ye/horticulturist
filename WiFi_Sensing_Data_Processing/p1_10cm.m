clear;
clc;


dBm_rate = xlsread('F:\Summer-Internship\raw data\p1_10cm.xlsx');
power_rate  = 10.^(dBm_rate/10-3);

repeat_num = 20;
test_num = 6;



delta_based_phase = zeros(repeat_num,test_num);
ratio_based_phase = zeros(repeat_num,test_num);
delta_based_amplitude = zeros(repeat_num,test_num);
ratio_based_amplitude = zeros(repeat_num,test_num);

for test_number = 1: test_num
    for repeat_number = 1: repeat_num
        path_ahead = 'F:\Summer-Internship\raw data\10cm\';
        path_ahead_contrast = 'F:\Summer-Internship\raw data\10cm\00';
        path_between = '_';
        path_behind = '.dat';

        %实验组path
        if test_number<10 && repeat_number < 10
            path_exp = [path_ahead,'0',num2str(test_number),path_between,'0',num2str(repeat_number),path_behind];
        end
        if test_number<10 && repeat_number >= 10
            path_exp = [path_ahead,'0',num2str(test_number),path_between,num2str(repeat_number),path_behind];
        end
        if test_number>=10 && repeat_number < 10
            path_exp = [path_ahead,num2str(test_number),path_between, '0',num2str(repeat_number),path_behind];
        end
        if test_number>=10 && repeat_number >= 10
            path_exp = [path_ahead,num2str(test_number),path_between, num2str(repeat_number),path_behind];
        end
    

        %对照组path
        if repeat_number < 10
            path_con = [path_ahead_contrast,path_between,'0',num2str(repeat_number),path_behind];
        end
        if repeat_number >= 10
            path_con = [path_ahead_contrast,path_between,num2str(repeat_number),path_behind];
        end

        path_con = 'F:\Summer-Internship\raw data\10cm\00_01.dat';

        %文件读入部分，确定文件长度并读入文件
        csi_trace = read_bf_file(path_exp);
        csi_trace_air = read_bf_file(path_con);
        L = size(csi_trace,1);
        L_air = size(csi_trace_air,1);

        if L>=L_air
            length = L_air;
        else
            length = L;
        end

        %amplitude代表波幅，phase代表相位，这里进行初始化，将向量全部置零
        %Initialization
        amplitudeA = zeros(length,30);
        amplitudeB = zeros(length,30);
        amplitudeC = zeros(length,30);
        amplitudeA_air = zeros(length,30);
        amplitudeB_air = zeros(length,30);
        amplitudeC_air = zeros(length,30);
        phaseA = zeros(L,30);
        phaseB = zeros(L,30);
        phaseC = zeros(L,30);
        phaseA_air = zeros(L,30);
        phaseB_air = zeros(L,30);
        phaseC_air = zeros(L,30);
        delta_phase = zeros(L,30);

        %首先从csi信心中获取对应30个载波波段的csi信息中的相位与幅值信息
        %Extract amplitude and phase from csi


        for m = 1:length
            csi_entry = csi_trace{m};
            csi = get_scaled_csi(csi_entry);
            csi_entry_air = csi_trace_air{m};
            csi_air = get_scaled_csi(csi_entry_air);

            %csi1=db(abs(squeeze(csi)));  %power

            csi_amp = abs(squeeze(csi));       %amplitude
            amplitudeA(m,:) = csi_amp(1,1,:);
            amplitudeB(m,:) = csi_amp(1,2,:);
            %amplitudeC(m,:) = csi_amp(1,3,:);
            csi_amp_air = abs(squeeze(csi_air));       %air amplitude
            amplitudeA_air(m,:) = csi_amp_air(1,1,:);
            amplitudeB_air(m,:) = csi_amp_air(1,2,:);
            %amplitudeC_air(m,:) = csi_amp_air(1,3,:);
            phaseA = squeeze(csi(1,1,:));   %phase
            phaseB = squeeze(csi(1,2,:));
            %phaseC = squeeze(csi(1,3,:));
            phaseA_air = squeeze(csi_air(1,1,:));   %air phase
            phaseB_air = squeeze(csi_air(1,2,:));
            %phaseC_air = squeeze(csi_air(1,3,:));
            delta = phaseB-phaseA;        %phase difference
            delta_phase(m,:)=abs(angle(delta));
            delta_air = phaseB_air-phaseA_air;    %air phase difference
            delta_phase_discuss(m,:)=angle(phaseA-phaseA_air);
            angle_melon(m,:) = angle(phaseA);
            delta_phase_air(m,:)=abs(angle(delta_air));
            phi_q(m,:) = abs(delta_phase(m,:)-delta_phase_air(m,:));
        end

        %对波幅进行降噪

        %Amplitude Denoising
        fs=100; %Sampling rate
        wp=1.5; %ͨPass-band frequency
        ws=4;   %Stop-band frequency
        rp=1;   %Max fading on pass-band
        as=80;  %Min fading on stop-band

        %巴特沃斯滤波
        %Butterworth filtering
        Y = amplitudeA;
        Y2 = amplitudeB;
        YfreqDomain = fft(Y);
        YfreqDomain2 = fft(Y2);
        YfreqDomain(2:length,:)=0;
        YfreqDomain2(2:length,:)=0;
        amplitudeA = ifft(YfreqDomain);
        amplitudeB = ifft(YfreqDomain2);
        amplitude_processed = abs(wifi_butterworth(amplitudeA,fs,wp,ws,rp,as));
        amplitude_processed2 = abs(wifi_butterworth(amplitudeB,fs,wp,ws,rp,as));
        amplitude_ratio = amplitude_processed./amplitude_processed2;
        amplitude_processed_air1 = abs(wifi_butterworth(amplitudeA_air,fs,wp,ws,rp,as));
        amplitude_processed_air2 = abs(wifi_butterworth(amplitudeB_air,fs,wp,ws,rp,as));
        amplitude_ratio_air = amplitude_processed_air1./amplitude_processed_air2;
        %amplitude_quotient
        A_q = amplitude_ratio./amplitude_ratio_air;
        %plot(A_q(:,1),'r');
        %hold on;
        %mean函数求复数相位角，mean函数为求平均数，std函数求标准偏差
        delta_amp = amplitude_processed_air1-amplitude_processed;
        aratio = zeros (30,1);
        %avg1为air的各channel幅值平均，avg2为带物品的各channel幅值平均
        avg1 = mean(mean(amplitude_processed_air1(:,:)));
        avg2 = mean(mean(amplitude_processed(:,:)));
        avg3 = angle(phaseA);
        avg4 = angle(phaseA_air);





        %Phase denosiing---------------------------------------------
        phi_q(L:L+300,:)=mean(mean(phi_q));
        n=300;k=0;m=0;

        %Sliding window average
        for i=1:L
            m = m+1;
            for l=1:30
                for j=i:n+i-1
                    k = k+1;
                    W(k,l)=phi_q(j,l);
                end
                phi_q_p(m,l)=mean(W(:,l));
            end
            k=0;
        end
        phi_q_p(1,:)=phi_q_p(2,:);
        %最后得到的两个原始信息应为phi_q_p和A_q

        %artio与phase_std用于channel的选择
        %aratio记录了所有的相对于空气的对照组的实验组数据对比，标准方差的大小
        for i =1:30
            aratio(i,1) = std (A_q(:,i));
        end

        %通过对phi进行方差得到相关相位方差数据
        for k =1:30
            phase_std(k,1) = std (phi_q_p(:,k));
        end


        aratio=aratio';
        p=phase_std';
        m1=zeros(1,4);
        a_min=zeros(length,30);
        phase_min=zeros(length,30);
        m2=zeros(1,4);
        amp_min=zeros(length,30);
        p_min=zeros(length,30);

        %根据相位信息进行载波选取(选取方差最小的4个数据组)
        %Phase related subcarrier selection
        %进行4次如下操作，首先找到当前相位方差的最小值对应的下标，将对应的值数组取出，然后把最小值覆盖为最大值
        %这样重复，最终得到方差最小的4组数据

        m2(1)=find(p==min(p));%Find the index of the first minimum value
        p_min(:,1)=phi_q_p(1:length,m2(1));%Find this value accoding to its index
        p(m2(1))=max(p);%Filter this index by assigning the maximum value to it

        m2(2)=find(p==min(p));
        p_min(:,2)=phi_q_p(1:length,m2(2));
        p(m2(2))=max(p);

        m2(3)=find(p==min(p));
        p_min(:,3)=phi_q_p(1:length,m2(3));
        p(m2(3))=max(p);

        m2(4)=find(p==min(p));
        p_min(:,4)=phi_q_p(1:length,m2(4));
        p(m2(4))=max(p);

        %根据我们得到的phase下标取出相应的amplitude原始数据
        for q=1:4
            amp_min(:,q)=A_q(1:length,m2(q));
        end


        %Amplitude related subcarrier selection
        m1(1)=find(aratio==min(aratio), 1);%Find the index of the first minimum value
        a_min(:,1)=A_q(1:length,m1(1));%Find this value accoding to its index
        aratio(m1(1))=max(aratio);%Filter this index by assigning the maximum value to it

        m1(2)=find(aratio==min(aratio), 1);
        a_min(:,2)=A_q(1:length,m1(2));
        aratio(m1(2))=max(aratio);

        m1(3)=find(aratio==min(aratio), 1);
        a_min(:,3)=A_q(1:length,m1(3));
        aratio(m1(3))=max(aratio);

        m1(4)=find(aratio==min(aratio), 1);
        a_min(:,4)=A_q(1:length,m1(4));
        aratio(m1(4))=max(aratio);

        %根据我们得到的amplitude下标取出相应的phase
        for q=1:4
            phase_min(:,q)=phi_q_p(1:length,m1(q));
        end

        %进行最后的平均值求取
        sum1 = 0;
        sum2 = 0;
        sum3 = 0;
        sum4 = 0;
        for l=1:4
            p_min(1,l)=mean(p_min(:,l));
            amp_min(1,l)=mean(amp_min(:,l));
            phase_min(1,l)=mean(phase_min(:,l));
            a_min(1,l)=mean(a_min(:,l));
            sum1 = sum1 + p_min(1,l);
            sum2 = sum2 + amp_min(1,l);
            sum3 = sum3 + phase_min(1,l);
            sum4 = sum4 + a_min(1,l);
        end

        %最后得到平均波幅与相位信息
        avg_delta_phase_based_phase = sum1/4;
        avg_amp_ratio_based_phase = sum2/4;
        avg_delta_phase_based_amplitude = sum3/4;
        avg_amp_ratio_based_amplitude = sum4/4;

        %power_caused_quotient = sqrt(power_rate(repeat_number,test_number+1)/power_rate(repeat_number,1));

        delta_based_phase(repeat_number,test_number)=avg_delta_phase_based_phase;
        ratio_based_phase(repeat_number,test_number)=avg_amp_ratio_based_phase;
        %ratio_based_phase(repeat_number,test_number)=avg_amp_ratio_based_phase/power_caused_quotient;
        delta_based_amplitude(repeat_number,test_number)=avg_delta_phase_based_amplitude;
        ratio_based_amplitude(repeat_number,test_number)=avg_amp_ratio_based_amplitude;
        %ratio_based_amplitude(repeat_number,test_number)=avg_amp_ratio_based_amplitude/power_caused_quotient;

    end
end

delta_based_phase = sort(delta_based_phase);
ratio_based_phase = sort(ratio_based_phase);
delta_based_amplitude = sort(delta_based_amplitude);
ratio_based_amplitude = sort(ratio_based_amplitude);

figure(11);
subplot(2,2,1);
plot(ratio_based_phase(:,1),'Linewidth',2);
hold on
plot(ratio_based_phase(:,2),'Linewidth',2);
hold on
plot(ratio_based_phase(:,3),'Linewidth',2);
hold on
plot(ratio_based_phase(:,4),'Linewidth',2);
hold on
plot(ratio_based_phase(:,5),'Linewidth',2);
hold on
plot(ratio_based_phase(:,6),'Linewidth',2);
hold on
plot(delta_based_phase(:,7),'Linewidth',2);
hold on
xlabel('实验编号');
legend('满水','空瓶');
title('1号天线测量幅值数据');

subplot(2,2,2);
plot(ratio_based_amplitude(:,1),'Linewidth',2);
hold on
plot(ratio_based_amplitude(:,2),'Linewidth',2);
hold on
plot(ratio_based_amplitude(:,3),'Linewidth',2);
hold on
plot(ratio_based_amplitude(:,4),'Linewidth',2);
hold on
plot(ratio_based_amplitude(:,5),'Linewidth',2);
hold on
plot(ratio_based_amplitude(:,6),'Linewidth',2);
hold on
plot(ratio_based_amplitude(:,7),'Linewidth',2);
hold on
xlabel('实验编号');
legend('满水','空瓶');
title('1号天线测量幅值数据');


%进行绘图


% beta_air = 2*pi*5*10^9 / 3*10^8
% alpha_air = 0000000000

%Compared amplitude methods----------------------------------------
%FFT denoising. [Unless I use apparant amplitude, this denoising is not selected.]
% Y = amplitude_processed(50:L,:);
% Y2 = amplitudeA(50:L,:);
% YfreqDomain = fft(Y);
% YfreqDomain2 = fft(Y2);
% % stem(abs(YfreqDomain)); %use abs command to get the magnitude
% % %similary, we would use angle command to get the phase plot!
% % %we'll discuss phase in another post though!
% % xlabel('Sample Number')
% % ylabel('Amplitude')
% % title('Using the Matlab fft command')
% % grid
% % axis([0,100,0,120])
% YfreqDomain(2:L,:)=0;
% YfreqDomain2(2:L,:)=0;
% Y = ifft(YfreqDomain);
% Y2 = ifft(YfreqDomain2);
% hold on
% plot(Y(:,1));
% title('FFT on original&butter'); %(similar with FFT on butterworth signal, difference is very slight)
%
% %DWT denoising[Comparison]
% [c,s] = wavedec2(amplitudeA,2,'coif3');
% [c2,s2] = wavedec2(amplitude_processed,2,'coif3');
% n=[1,2];
% p=[10.12,23.28];
% nc = wthcoef2('h',c,s,n,p,'s');
% nc = wthcoef2('v',nc,s,n,p,'s');
% nc = wthcoef2('d',nc,s,n,p,'s');
% x1 = waverec2(nc,s,'coif3');
% nc2 = wthcoef2('h',c2,s2,n,p,'s');
% nc2 = wthcoef2('v',nc2,s2,n,p,'s');
% nc2 = wthcoef2('d',nc2,s2,n,p,'s');
% x2 = waverec2(nc2,s,'coif3');
% subplot(5,1,4);
% %plot(ca1); title('low frequency info of 1st wavedec');
% plot(x1(:,1)); title('DWT reconstructed signal(on original)');
% subplot(5,1,5);
% plot(x2(:,1)); title('DWT reconstructed signal(on butterworth)');
% print('Denoising comparison','-dpng');

%Butterworth processed VS. Amplitude ratio, std plot
% hold on
% plot(aratio,'r','Linewidth',2.5);
% hold on
% plot(bratio,'b','Linewidth',2.5);
% hold on
% plot(cratio,'g','Linewidth',2.5);
% xlabel('Subcarrier Index');
% ylabel('Amplitude Variance');
% L=legend('Neighboring antennas','Antenna 1','Antenna 2');
% set(L,'FontName','Times New Roman','FontSize',11,'FontWeight','normal');
%print('Amplitude Ratio','-dpng');
% %Butterworth processed result visualization
% plot(amplitude_ratio(:,1)); title('Ratio denoised');