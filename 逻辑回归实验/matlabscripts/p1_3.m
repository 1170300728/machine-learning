clear;
% 
% %构造一个初始边界theta1，对应theta，并画出该线
% theta1=[20,10,-130];%这里随便写了一个初始边界
% x1=0:0.001:10;
% y1=(theta1(2).*x1+theta1(3))./(-theta1(1));
% hold on
% plot(x1,y1)
% 
% %选取10*10的样本点，分别位于（1,1）到（10,10），并加入高斯噪声成为x数据集
% m=100;
% x=zeros(3,100);
% x(3,:)=ones(1,100);
% for i=1:10
%     for j=1:10
%         x(1,(i-1)*10+j)=j;
%         x(2,(i-1)*10+j)=i;
%     end
% end
% z1=0.5*randn(1,m); %产生高斯噪声
% x(1,:)=x(1,:)+z1;
% z2=0.5*randn(1,m); %产生高斯噪声
% x(2,:)=x(2,:)+z2;
% 
% %构造y数据集，利用初始边界theta1，判断x中各点位置处于边界哪一册并确定y值
% y=zeros(100,1);
% for i=1:100
%     if(x(1,i)*theta1(1)+x(2,i)*theta1(2)+theta1(3)>0)
%         y(i)=1;
%         scatter(x(2,i),x(1,i),'r');
%     else
%         scatter(x(2,i),x(1,i),'b');
%     end
% end

x=load("haberman.data");
y=x(:,4)-ones(size(x,1),1);
m=size(x,1);
for i=1:m
    if(x(i,4)==2)
%        scatter3(x(i,1),x(i,2),x(i,3),"b");
        x(i,4)=1;
    else
    end
end
y1=sum(y(:)==0);
x1=zeros(y1,4);
x2=zeros(numel(y)-y1,4);
i1=1;i2=1;
for i=1:m
    if(y(i)==0)
        x1(i1,:)=x(i,:);
        i1=i1+1;
    else
        x2(i2,:)=x(i,:);
        i2=i2+1;
    end
end
    subplot(1,2,1)%将窗口两等分，并作图
        scatter3(x1(:,1),x1(:,2),x1(:,3),"r");
        hold on
        scatter3(x2(:,1),x2(:,2),x2(:,3),"b");
        hold on
xx=x';
x=zeros(3,306);
x(1,:)=xx(1,:);
x(2,:)=xx(3,:);
x(3,:)=xx(4,:);
    subplot(1,2,2)%将窗口两等分，并作图
        scatter(x1(:,1),x1(:,3),"r");
        hold on
        scatter(x2(:,1),x2(:,3),"b");
        hold on
%梯度下降法
theta=[1;1;1];%初始化theta
alpha=0.001;%步长
temptheta=[100,100,100];
k=0;%计迭代次数
lamda=0.001;%正则项系数
l=1000000;
delta=zeros(1,l);
while abs(sum(temptheta-theta))>10^-6 %计算这次与上次的theta值的差，若足够小则认为已经完成逻辑回归
% while k<l
%     temptheta=theta;
%     %根据推导的公式计算迭代过程
%     jd= x*(1./(1+exp(-x'*temptheta))-y)./100;
%     h=zeros(1,m);
%     for i=1:m
%         h(i)=1./(1+exp(-x(:,i)'*temptheta));
%     end
%     
%     jdd=x*diag(h)*x'./m;
%     
%     theta=temptheta-(inv(jdd)*jd+lamda*theta).*alpha;
%     k=k+1;
    
    temptheta=theta;
    %根据推导的公式计算迭代过程
    %jd= x*(1./(1+exp(-x'*temptheta))-y)./m;%无正则项
    jd= x*(1./(1+exp(-x'*temptheta))-y)./m+lamda.*temptheta;%有正则项
    theta=temptheta-jd.*alpha;
    k=k+1;
    delta(k)=abs(sum(temptheta-theta));
end

fprintf("k=%d\ndeltatheta=%d\n",k,abs(sum(temptheta-theta)));

%根据theta画出最终的逻辑回归边界
y0=0:1:100;
x0=((theta(1).*y0)+theta(3))./(-theta(2));
plot(y0,x0,"bla");
suptitle('梯度下降法实现逻辑回归')