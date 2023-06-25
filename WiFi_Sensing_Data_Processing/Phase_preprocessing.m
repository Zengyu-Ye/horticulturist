clear;
clc;

%Phase denoising
csi_trace=read_bf_file('./csi/isabelle/1223/melon3.dat');
%csi_trace2=read_bf_file('./csi/isabelle/experiment/apple_2.dat');
L= length(csi_trace);

phaseA=zeros(L,30);
phaseB=zeros(L,30);
phaseC=zeros(L,30);
delta_phase=zeros(L,30);

for m=1:L
    csi_entry=csi_trace{m};
    %csi_entry2=csi_trace2{m};
    csi=get_scaled_csi(csi_entry);
    %csi2=get_scaled_csi(csi_entry2);
    phaseA = squeeze(csi(1,1,:));
    %phaseA2 = squeeze(csi(2,1,:));
    phaseB = squeeze(csi(1,2,:));
    phaseC = squeeze(csi(1,3,:));
    delta = phaseC-phaseA;
    delta_phase(m,:)=abs(angle(delta));
    phase_origin(m,:)=abs(angle(phaseA));
end

delta_phase(L:L+300,:)=mean(mean(delta_phase));   
phase_origin(900,30)=0; 
n=300;k=0;m=0;

for i=1:L
    m = m+1;
    for l=1:30
        for j=i:n+i-1
            k = k+1;
            W(k,l)=delta_phase(j,l);
        end
    delta_phase_p(m,l)=mean(W(:,l));      
    end
k=0;
end
delta_phase_p(1,:)=delta_phase_p(2,:);


%     phase_origin_p(j,:)=mean(phase_origin(300*(j-1)+1:300*j,:));
% end

%Subcarrier selection
for k =1:30
    phase_std(k,1) = std (delta_phase_p(:,k));
end

a=phase_std';
m=zeros(1,4);                                      
a_min=zeros(L,30); 
m(1)=find(a==min(a));%Find the index of the first minimum value
a_min(:,1)=delta_phase_p(:,m(1));%Find this value accoding to its index
a(m(1))=max(a);%Filter this index by assigning the maximum value to it

m(2)=find(a==min(a));
a_min(:,2)=delta_phase_p(:,m(2));
a(m(2))=max(a);

m(3)=find(a==min(a));
a_min(:,3)=delta_phase_p(:,m(3));
a(m(3))=max(a);

m(4)=find(a==min(a));
a_min(:,4)=delta_phase_p(:,m(4));
a(m(4))=max(a);

sum = 0;
for l=1:4
    a_min(1,l)=mean(a_min(:,l));
    sum = sum + a_min(1,l);
end
avg_delta_phase = sum/4;



% hold on;
% polar(phase_origin,'b');
% polar(delta_phase,'r');
% % plot(real(phaseA),imag(phaseA),'r');
% % plot(real(delta),imag(delta),'b');
% legend('Original phases','Processed phases between antenna1-2');
plot(delta_phase_p(:,1));
% plot(phase_std,'m');
% hold on;
% plot(phase_std_origin,'y');
% xlabel('Subcarrier Index');
% ylabel('Phase Difference Variance');
% % print('Subcarrier selection','-dpng');
