clear all;
n=1;%特征向量个数
m=20;%样本个数
%建立x和y矩阵，用于保存函数。
x=ones(n,m);
y=ones(n,m);

for i=1:n
    x(i,:)=linspace(0,10,m);%从0到10的区间平均采样
    y(i,:)=sin(pi/2*x(i,:));%计算样本值
    z1=0.35*randn(1,m); %产生高斯噪声，取较小系数减少噪声否则过于难拟合
    y(i,:)=y(i,:)+z1; %将高斯噪声加入正弦波中
    scatter(x(i,:),y(i,:))%画图
    hold on
end

[~,k]=size(x);%获取x的列数
i=1;
X0=zeros(i+1,k);%将矩阵初始化为0
for k0=1:k
    for n0=1:i+1
        X0(n0,k0)=x(k0)^(n+1-n0);%循环构造矩阵X0
    end
end
X=X0';%按推导要求对矩阵X0求转置
w=[0;0];%先简化问题为一阶线性函数拟合，需要常数项和一次项
theta=0.0001;%取步长
eta=0.001;%取系数
for j=1:50%设置最大训练次数
    w0=X'*X*w-X'*y'+theta*w;
    w=w-eta*w0;
end

%这一部分为了画出更平滑的曲线取得x值相当密集
x0=0:0.001:10;
y0=x0*w(1,1)+w(2,1);
plot(x,y,'*')
hold on
plot(x0,y0)

suptitle('一阶线性函数拟合正弦曲线')