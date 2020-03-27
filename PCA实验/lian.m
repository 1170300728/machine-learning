clear;

png='.png';
jpg='.jpg';
k=1;
data=zeros(60,900);
for i=396:455
    file=strcat(num2str(i),jpg);
     data30=imread(file);
     if k<=10
     subplot(2,10,k);
     imshow(data30);
     hold on
     end
    data30=double(rgb2gray(data30));
    data(k,:)=reshape(data30,[900,1]);
    k=k+1;
end

[dotNum,dimention]=size(data);
k=20;


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


Pn=sum(sum(data.^2));
Pk=sum(sum((data-kData+meanData).^2));
snr=10*log10(Pn/Pk)
for i=1:10
    kData30=reshape(kData(i,:),[30,30]);
    kData30=uint8(kData30);
    subplot(2,10,i+10);
    imshow(kData30);
    hold on
end