clear;

%构造一个初始边界theta1，对应theta，并画出该线
theta1=[20,10,-130];%这里随便写了一个初始边界
x1=0:0.001:10;
y1=(theta1(2).*x1+theta1(3))./(-theta1(1));
hold on
plot(x1,y1)

%选取10*10的样本点，分别位于（1,1）到（10,10），并加入高斯噪声成为x数据集
m=100;
x=zeros(3,100);
x(3,:)=ones(1,100);
for i=1:10
    for j=1:10
        x(1,(i-1)*10+j)=j;
        x(2,(i-1)*10+j)=i;
    end
end
z1=0.5*randn(1,m); %产生高斯噪声
x(1,:)=x(1,:)+z1;
z2=0.5*randn(1,m); %产生高斯噪声
x(2,:)=x(2,:)+z2;

%构造y数据集，利用初始边界theta1，判断x中各点位置处于边界哪一册并确定y值
y=zeros(100,1);
for i=1:100
    if(x(1,i)*theta1(1)+x(2,i)*theta1(2)+theta1(3)>0)
        y(i)=1;
        scatter(x(2,i),x(1,i),'r');
    else
        scatter(x(2,i),x(1,i),'b');
    end
end

%梯度下降法
theta=[0;0;0];%初始化theta
alpha=0.1;%步长
temptheta=[100,100,100];
k=0;%计迭代次数
lamda=0.001;
while abs(sum(temptheta-theta))>10^-15 %计算这次与上次的theta值的差，若足够小则认为已经完成逻辑回归
    temptheta=theta;
    %根据推导的公式计算迭代过程
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

%根据theta画出最终的逻辑回归边界
x0=0:0.001:10;
y0=(theta(2).*x0+theta(3))./(-theta(1));
hold on
plot(x0,y0,"g")

suptitle('牛顿法实现逻辑回归')