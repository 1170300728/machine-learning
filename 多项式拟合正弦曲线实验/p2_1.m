clear all;
n=1;%������������
m=20;%��������
%����x��y�������ڱ��溯����
x=ones(n,m);
y=ones(n,m);
theta=zeros(n,1);%���ò���,n��1�е�������������theta

for i=1:n
    x(i,:)=linspace(0,10,m);%��0��10������ƽ������
    y(i,:)=sin(pi/2*x(i,:));%��������ֵ
    z1=0.35*randn(1,m); %������˹������ȡ��Сϵ����������������������
    y(i,:)=y(i,:)+z1; %����˹�����������Ҳ���
    scatter(x(i,:),y(i,:))%��ͼ
    hold on
end

[~,k]=size(x);
for i=1:9
    X0=zeros(i+1,k);%�������ʼ��Ϊ0
    for k0=1:k
        for n0=1:i+1
            X0(n0,k0)=x(k0)^(i+1-n0);%ѭ���������X0
        end
    end
    X=X0';%���Ƶ�Ҫ��Ծ���X0��ת��
    ANSS=(X'*X)\X'*y';%���ݹ�ʽ���
    for j=1:i+1
        theta(j,i)=ANSS(j);%answer����洢ÿ����õķ���ϵ�������д洢
    end
    x0=0:0.001:10;
    y0=ANSS(1)*x0.^i    ;%������õ�ϵ����ʼ�����������ʽ����
    for num=2:1:i+1
        y0=y0+ANSS(num)*x0.^(i+1-num);%����׼����ͼ
    end
    subplot(3,3,i)%�����ھŵȷ֣�����ͼ
    plot(x,y,'*')
    hold on
    plot(x0,y0)
end
suptitle('��ͬ��������������Ͻ������1��9��')