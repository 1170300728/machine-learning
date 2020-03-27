clear all;
n=1;%特征向量个数
m=10;%样本个数
%建立x和y矩阵，用于保存函数。
x=ones(n,m);
y=ones(n,m);

for i=1:n
    x(i,:)=linspace(0,4,m);%从0到10的区间平均采样
    y(i,:)=sin(pi/2*x(i,:));%计算样本值
    z1=0.05*randn(1,m); %产生高斯噪声，取较小系数减少噪声否则过于难拟合
    y(i,:)=y(i,:)+z1; %将高斯噪声加入正弦波中
    scatter(x(i,:),y(i,:))%画图
    hold on
end

[~,k]=size(x);%获取x的列数
i=3;
X0=zeros(i+1,k);%将矩阵初始化为0
for k0=1:k
    for n0=1:i+1
        X0(n0,k0)=x(k0)^(i+1-n0);%构造矩阵X0
    end
end
X=X0';%按推导要求对矩阵X0求转置
w=zeros(i+1,1);%先简化问题为一阶线性函数拟合，需要常数项和一次项
theta=0.01;%取步长
eta=0.00005;%取系数
for j=1:1000000%设置最大训练次数
    w0=X'*X*w-X'*y'+theta*w;
    w=w-eta*w0;
end

%这一部分为了画出更平滑的曲线取得x值相当密集
x0=0:0.001:4;
y0=zeros(1,size(x0,2));
for j=0:i
    y0=x0.^j*w(i+1-j,1)+y0;
end
hold on
plot(x0,y0)

suptitle('梯度下降法高阶线性函数拟合曲线');