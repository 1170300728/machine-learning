clear;
%��ʼ��ģ�Ͳ���
dimention=3;
k=2;
dotNum=100;

data=zeros(dotNum,dimention);
data(:,1)=repmat(1:1:10,1,10);
for i=1:dotNum
    data(i,2)=ceil(i/10);
end
data(:,3)=data(:,1)+data(:,2);
data=data+randn(dotNum,dimention);
scatter3(data(:,1),data(:,2),data(:,3),'g');
hold on

meanData=ones(size(data,1),1)*mean(data);
data=data-meanData;
dataCov=cov(data);

[V,D]=eig(dataCov);
for i=1:dimention
    for j=i+1:dimention
        if D(i,i)<D(j,j)
            tempD=D(i,i);
            D(i,i)=D(j,j);
            D(j,j)=tempD;
            tempV=V(:,i);
            V(:,i)=V(:,j);
            V(:,j)=tempV;
        end
    end
end
kD=D(:,1:k);
kV=V(:,1:k);
kData=data*kV*kV'+meanData;
scatter3(kData(:,1),kData(:,2),kData(:,3),'r');
view(45,30);
% %ÿ����һ������
% %newX  ��ά����¾���
% %T �任����
% %meanValue  Xÿ�о�ֵ���ɵľ������ڽ���ά��ľ���newX�ָ���X
% %CRate ������
% %�������Ļ���������
% meanValue=ones(size(X,1),1)*mean(X);
% X=X-meanValue;%ÿ��ά�ȼ�ȥ��ά�ȵľ�ֵ
% C=cov(x);%����Э�������
%  
% %������������V������ֵD
% [V,D]=eig(C);
% %��������������������
% [dummy,order]=sort(diag(D),'descend');
% V=V(:,order);%������������������ֵ��С���н�������
% d=diag(D);%������ֵȡ��������һ��������
% newd=d(order);%������ֵ���ɵ�����������������
%  
% %ȡǰn���������������ɱ任����
% sumd=sum(newd);%����ֵ֮��
% for j=1:length(newd)
%     i=sum(newd(1:j,1))/sumd;%���㹱���ʣ�������=ǰn������ֵ֮��/������ֵ֮��
%     if i>CRate%�������ʴ���95%ʱѭ������,������ȡ���ٸ�����ֵ
%         cols=j;
%         break;
%    end
% end
% T=V(:,1:cols);%ȡǰcols���������������ɱ任����T
% newX=X*T;%�ñ任����T��X���н�ά
