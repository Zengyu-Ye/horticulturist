amplitude=amplitudeA;
[r,c]=size(amplitude);%��ȡ��������������
me=mean(amplitude);%mean(A)�����еĸ�����Ϊ����
m=repmat(me,[r,1]);%��me�����ݶѵ��ڣ�rx1���ľ���m��
B=amplitude-m;%ÿ�����ȥƽ��ֵ
C=cov(B);%���ش�СΪM*N����cov(X)��СΪN*N�ľ���cov(X)�ĵ�(i,j)��Ԫ�ص���X�ĵ�i���������j�������ķ����C(Xi,Xj)
[ev,ed]=eig(C);%�����A��ȫ������ֵ�����ɶԽ���ed������A��������������ev����������

M=ev(:,30);
out=amplitude*M;
fangcha=var(out);
%plot(out,'r');
tt1 = out(60:170,:);
tt2 = out(360:430,:);
%hold on

M=ev(:,29);
out=amplitude*M;
%out=diff(out);
out = abs(out);
%plot(out,'b')
fangcha1=var(out);

 M=ev(:,28);
 out1=amplitude*M;
 fangcha2=var(out1);
% plot(out1,'g')
% 
% M=ev(:,27);
% out=amplitudeB*M;
% plot(out,'y')
% %legend('first','second','third','fourth')
% legend('first','second','third')
% title('PCA')
% hold off