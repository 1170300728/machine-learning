clear;
%��ʼ�����ģ������
gaussNum=4;%��˹ģ�͸���
dimention=2;%ģ��ά��
dotNum=1000;%�����㼯��С
width=15;%����ȡֵ��Χ

%��˹��ϳ�ʼ��
centers=width*rand(dimention,gaussNum);%��˹ģ������
weight = rand(1,gaussNum);%ģ�Ͳ�����������
weight = weight / norm(weight, 1); % ��һ������֤������Ϊ1
weight = (weight+0.2)./(gaussNum*0.2+1);%��������ʹ�������ٵĲ��ֲ�����
%�����ģ�Ͱ����������ʵ������
n=zeros(1,gaussNum);%���������
for i = 1: gaussNum
    if i ~= gaussNum
        n(i) = floor(dotNum*weight(i));
    else
        n(i) = dotNum - sum(n);
    end
end
%�������˹ģ�͵�Э�������
Cs=zeros(dimention,dimention,gaussNum);%Э�������
for i=1:gaussNum
    for j=1:dimention
        Cs(j,j,i)=sqrt((rand(1)+0.2)*2);%Э�������Խ���
        for k=1:dimention
            if j~=k
                Cs(j,k,i)=(floor(rand(1)*8)-4);%Э�������ǶԽ���Ԫ��
                if Cs(j,k,i)<0%������������������
                    Cs(j,k,i)=-sqrt(10-Cs(j,k,i))/10;
                else
                    Cs(j,k,i)=sqrt(10+Cs(j,k,i))/10;
                end
                Cs(k,j,i)=Cs(j,k,i);%Э���������һ���Գ��󣬶Խ�λ����ֵ��ͬ
            end
        end
    end
end

%���
data=zeros(dimention,dotNum);%ʵ�ʲ����㼯
start = 1;
for i=1:gaussNum
    data(:,start:start+n(i)-1)=repmat(centers(:,i),1,n(i))+Cs(:,:,i)*randn(dimention,n(i));
    start=start+n(i);
end

%��ͼ������ɫȦ�������в����㣬���ú�ɫ*���ģ��ʵ������
for i=1:gaussNum
    scatter(centers(1,i),centers(2,i),"*",'bla');
    hold on
end
for i=1:dotNum
    scatter(data(1,i),data(2,i),'g');
    hold on
end
%%
%k-means�������࣬��ý�Ϊ׼ȷ�ĳ�ʼλ��
%��ʼ����������
maxStepsm=200;%����������
kMeans=width*rand(dimention,gaussNum);%������ɾ����ʼ����
kMeans0=kMeans;%��¼��һ�εľ�������
kMeansRecord=zeros(dimention*(maxStepsm+1),gaussNum);%��¼�������
%��¼��һ��λ��
recordLine=0;
kMeansRecord(recordLine+1:recordLine+dimention,:)=kMeans;
m=0;%��������
while m<maxStepsm
    %��ʼ����������
    kMeanCenters=zeros(dimention+1,gaussNum);
    distances=zeros(gaussNum,dotNum);%ÿ�������㵽�ϴε�����ÿ���������ĵľ���
    shortcuts=zeros(1,dotNum);%ÿ�������㵽����ľ������ĵľ���
    nearCenters=zeros(1,dotNum);%����ÿ������������ľ������ı��
    %�������
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
    %ͨ����������µľ�������
    for i=1:dotNum
        kMeanCenters(1:dimention,nearCenters(i))=(kMeanCenters(1:dimention,nearCenters(i))*kMeanCenters(1+dimention,nearCenters(i))+data(:,i))./(kMeanCenters(1+dimention,nearCenters(i))+1);
        kMeanCenters(1+dimention,nearCenters(i))=kMeanCenters(1+dimention,nearCenters(i))+1;
    end
    %��¼����Ǩ��
    kMeans=kMeanCenters(1:dimention,:);
    recordLine=recordLine+dimention;
    kMeansRecord(recordLine+1:recordLine+dimention,:)=kMeans;
    m=m+1;
    %��������д���ȫ��㣬˵������������и����Ĳ����κβ�������������������µ�λ�ù����µ����ļ������
    for i=1:gaussNum
        if kMeans(:,i)==zeros(2,1)
            kMeans(:,i)=width*rand(dimention,1);
        end
    end
    %������ε���ֵ��ȫһ�£��˳�ѭ��
    if kMeans==kMeans0
        break;
    end
    kMeans0=kMeans;
end
%%
%��ͼ
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
%EM�㷨�����׼ȷ�ľ�������
maxStepsn=200;
data=data';
EMs=kMeans';
EMs0=EMs;%��¼��һ�εľ�������
EMsRecord=zeros(dimention*(maxStepsn+1),gaussNum);%��¼�������
%��¼��һ��λ��
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
    
    prep = sum(p,1)./size(p,1); %%��p��ÿһ�м��������ܵõ�ÿһ��������������
    
    EMs = p'*data; %%�ֱ�õ�dataÿһά����ÿһ�������������EMs(i,j),i��ά����j�Ǿ�����
    EMs= EMs./repmat((sum(p,1))',1,size(EMs,2));
    
    for j = 1 : length(prep)
        vari = repmat(p(:,j),1,size(data,2)).*(data- repmat(EMs(j,:),size(data,1),1)); %%�õ�һ���ض�����Xÿһά�ķ�����󣬳���p�����൱��ѡ������ڸþ����d�����㣩
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
%��ͼ
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