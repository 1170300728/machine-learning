clear;
%初始化混合模型数据
gaussNum=4;%高斯模型个数
dimention=2;%模型维数
dotNum=1000;%采样点集大小
width=15;%中心取值范围

%高斯混合初始化
centers=width*rand(dimention,gaussNum);%高斯模型中心
weight = rand(1,gaussNum);%模型采样数量比例
weight = weight / norm(weight, 1); % 归一化，保证比例合为1
weight = (weight+0.2)./(gaussNum*0.2+1);%调整比例使比例最少的部分不过低
%构造各模型包含采样点的实际数量
n=zeros(1,gaussNum);%采样点个数
for i = 1: gaussNum
    if i ~= gaussNum
        n(i) = floor(dotNum*weight(i));
    else
        n(i) = dotNum - sum(n);
    end
end
%构造各高斯模型的协方差矩阵
Cs=zeros(dimention,dimention,gaussNum);%协方差矩阵
for i=1:gaussNum
    for j=1:dimention
        Cs(j,j,i)=sqrt((rand(1)+0.2)*2);%协方差矩阵对角线
        for k=1:dimention
            if j~=k
                Cs(j,k,i)=(floor(rand(1)*8)-4);%协方差矩阵非对角线元素
                if Cs(j,k,i)<0%根据正负做开方处理
                    Cs(j,k,i)=-sqrt(10-Cs(j,k,i))/10;
                else
                    Cs(j,k,i)=sqrt(10+Cs(j,k,i))/10;
                end
                Cs(k,j,i)=Cs(j,k,i);%协方差矩阵是一个对称阵，对角位置数值相同
            end
        end
    end
end

%混合
data=zeros(dimention,dotNum);%实际采样点集
start = 1;
for i=1:gaussNum
    data(:,start:start+n(i)-1)=repmat(centers(:,i),1,n(i))+Cs(:,:,i)*randn(dimention,n(i));
    start=start+n(i);
end

%作图，用绿色圈画出所有采样点，并用红色*标出模型实际中心
for i=1:gaussNum
    scatter(centers(1,i),centers(2,i),"*",'bla');
    hold on
end
for i=1:dotNum
    scatter(data(1,i),data(2,i),'g');
    hold on
end
%%
%k-means方法聚类，获得较为准确的初始位置
%初始化聚类数据
maxStepsm=200;%最大迭代次数
kMeans=width*rand(dimention,gaussNum);%随机生成聚类初始中心
kMeans0=kMeans;%记录上一次的聚类中心
kMeansRecord=zeros(dimention*(maxStepsm+1),gaussNum);%记录聚类过程
%记录第一个位置
recordLine=0;
kMeansRecord(recordLine+1:recordLine+dimention,:)=kMeans;
m=0;%迭代次数
while m<maxStepsm
    %初始化迭代数据
    kMeanCenters=zeros(dimention+1,gaussNum);
    distances=zeros(gaussNum,dotNum);%每个采样点到上次迭代中每个聚类中心的距离
    shortcuts=zeros(1,dotNum);%每个采样点到最近的聚类中心的距离
    nearCenters=zeros(1,dotNum);%距离每个采样点最近的聚类中心编号
    %计算距离
    for i=1:dotNum
        for j=1:gaussNum
            dis=0;
            for k=1:dimention
                dis=dis+(data(k,i)-kMeans(k,j))^2;
            end
            distances(j,i)=sqrt(dis);
        end
        shortcuts(i)=min(distances(:,i));
        for j=1:gaussNum
            if shortcuts(i)==distances(j,i)
                nearCenters(i)=j;
            end
        end
    end
    %通过距离计算新的聚类中心
    for i=1:dotNum
        kMeanCenters(1:dimention,nearCenters(i))=(kMeanCenters(1:dimention,nearCenters(i))*kMeanCenters(1+dimention,nearCenters(i))+data(:,i))./(kMeanCenters(1+dimention,nearCenters(i))+1);
        kMeanCenters(1+dimention,nearCenters(i))=kMeanCenters(1+dimention,nearCenters(i))+1;
    end
    %记录中心迁移
    kMeans=kMeanCenters(1:dimention,:);
    recordLine=recordLine+dimention;
    kMeansRecord(recordLine+1:recordLine+dimention,:)=kMeans;
    m=m+1;
    %如果中心中存在全零点，说明聚类过程中有个中心不与任何采样点距离近，重新随机新的位置构成新的中心加入迭代
    for i=1:gaussNum
        if kMeans(:,i)==zeros(2,1)
            kMeans(:,i)=width*rand(dimention,1);
        end
    end
    %如果两次迭代值完全一致，退出循环
    if kMeans==kMeans0
        break;
    end
    kMeans0=kMeans;
end
%%
%作图
m=m+1;
for i=1:gaussNum
    x=1:dimention:m*dimention-1;
    y=zeros(dimention,m);
    for j=1:dimention
        for k=1:m
            y(j,k)=kMeansRecord(x(k)+j-1,i);
        end
    end
    plot(y(1,:),y(2,:),'b');
    hold on
    scatter(kMeans(1,i),kMeans(2,i),'*','b');
    hold on
end
%%
%EM算法计算更准确的聚类中心
maxStepsn=200;
data=data';
EMs=kMeans';
EMs0=EMs;%记录上一次的聚类中心
EMsRecord=zeros(dimention*(maxStepsn+1),gaussNum);%记录聚类过程
%记录第一个位置
recordLine=0;
EMsRecord(recordLine+1:recordLine+dimention,:)=EMs';
cov=zeros(dimention,dimention,gaussNum);
for i=1:gaussNum
    cov(:,:,i) = eye(dimention);
end
p=zeros(dotNum,gaussNum);
prep=ones(gaussNum)./gaussNum;
n=0;
while n<maxStepsn
    for i=1:gaussNum
        p(:,i)=prep(i)*mvnpdf(data,EMs(i,:),cov(:,:,i));
    end
    p=p./repmat(sum(p,2),1,size(p,2));
    
    prep = sum(p,1)./size(p,1); %%把p的每一行加起来就能得到每一个聚类的先验概率
    
    EMs = p'*data; %%分别得到data每一维对于每一个聚类的期望，EMs(i,j),i是维数，j是聚类数
    EMs= EMs./repmat((sum(p,1))',1,size(EMs,2));
    
    for j = 1 : length(prep)
        vari = repmat(p(:,j),1,size(data,2)).*(data- repmat(EMs(j,:),size(data,1),1)); %%得到一个特定聚类X每一维的方差矩阵，乘以p，（相当于选择出属于该聚类的d采样点）
        cov(:,:,j) = (vari'*vari)/sum(p(:,j),1);
    end
    if n>5
        if sum(abs(EMs0-EMs))<10^-7
            break;
        end
    end
    EMs0=EMs;
    recordLine=recordLine+dimention;
    EMsRecord(recordLine+1:recordLine+dimention,:)=EMs';
    n=n+1;
end
%%
%作图
n=n+1;
for i=1:gaussNum
    x=1:dimention:n*dimention-1;
    y=zeros(dimention,n);
    for j=1:dimention
        for k=1:n
            y(j,k)=EMsRecord(x(k)+j-1,i);
        end
    end
    plot(y(1,:),y(2,:),'r');
    hold on
    scatter(EMs(i,1),EMs(i,2),'*','r');
    hold on
end