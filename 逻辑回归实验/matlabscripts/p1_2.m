clear;

%����һ����ʼ�߽�theta1����Ӧtheta������������
theta1=[20,10,-130];%�������д��һ����ʼ�߽�
x1=0:0.001:10;
y1=(theta1(2).*x1+theta1(3))./(-theta1(1));
hold on
plot(x1,y1)

%ѡȡ10*10�������㣬�ֱ�λ�ڣ�1,1������10,10�����������˹������Ϊx���ݼ�
m=100;
x=zeros(3,100);
x(3,:)=ones(1,100);
for i=1:10
    for j=1:10
        x(1,(i-1)*10+j)=j;
        x(2,(i-1)*10+j)=i;
    end
end
z1=0.5*randn(1,m); %������˹����
x(1,:)=x(1,:)+z1;
z2=0.5*randn(1,m); %������˹����
x(2,:)=x(2,:)+z2;

%����y���ݼ������ó�ʼ�߽�theta1���ж�x�и���λ�ô��ڱ߽���һ�Ტȷ��yֵ
y=zeros(100,1);
for i=1:100
    if(x(1,i)*theta1(1)+x(2,i)*theta1(2)+theta1(3)>0)
        y(i)=1;
        scatter(x(2,i),x(1,i),'r');
    else
        scatter(x(2,i),x(1,i),'b');
    end
end

%�ݶ��½���
theta=[0;0;0];%��ʼ��theta
alpha=0.1;%����
temptheta=[100,100,100];
k=0;%�Ƶ�������
lamda=0.001;
while abs(sum(temptheta-theta))>10^-15 %����������ϴε�thetaֵ�Ĳ���㹻С����Ϊ�Ѿ�����߼��ع�
    temptheta=theta;
    %�����Ƶ��Ĺ�ʽ�����������
    jd= x*(1./(1+exp(-x'*temptheta))-y)./100;
    h=zeros(1,m);
    for i=1:m
        h(i)=1./(1+exp(-x(:,i)'*temptheta));
    end
    
    jdd=x*(diag(h)+0.00001)*x'./m;
    
    theta=temptheta-(inv(jdd)*jd+lamda*theta).*alpha;
    k=k+1;
end
fprintf("k=%d\n",k);

%����theta�������յ��߼��ع�߽�
x0=0:0.001:10;
y0=(theta(2).*x0+theta(3))./(-theta(1));
hold on
plot(x0,y0,"g")

suptitle('ţ�ٷ�ʵ���߼��ع�')