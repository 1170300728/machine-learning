clear all;
n=1;%特征向量个数
m=20;%样本个数
%建立x和y矩阵，用于保存函数。
x=ones(n,m);
y=ones(n,m);
theta=zeros(n,1);%设置参数,n行1列的向量，保存结果theta

for i=1:n
    x(i,:)=linspace(0,10,m);%从0到10的区间平均采样
    y(i,:)=sin(pi/2*x(i,:));%计算样本值
    z1=0.35*randn(1,m); %产生高斯噪声，取较小系数减少噪声否则过于难拟合
    y(i,:)=y(i,:)+z1; %将高斯噪声加入正弦波中
    scatter(x(i,:),y(i,:))%画图
    hold on
end

[~,k]=size(x);
for i=1:9
    X0=zeros(i+1,k);%将矩阵初始化为0
    for k0=1:k
        for n0=1:i+1
            X0(n0,k0)=x(k0)^(i+1-n0);%循环构造矩阵X0
        end
    end
    X=X0';%按推导要求对矩阵X0求转置
    ANSS=(X'*X)\X'*y';%根据公式求解
    for j=1:i+1
        theta(j,i)=ANSS(j);%answer矩阵存储每次求得的方程系数，按列存储
    end
    x0=0:0.001:10;
    y0=ANSS(1)*x0.^i    ;%根据求得的系数初始化并构造多项式方程
    for num=2:1:i+1
        y0=y0+ANSS(num)*x0.^(i+1-num);%采样准备作图
    end
    subplot(3,3,i)%将窗口九等分，并作图
    plot(x,y,'*')
    hold on
    plot(x0,y0)
end
suptitle('不同次数方程曲线拟合结果，从1到9阶')