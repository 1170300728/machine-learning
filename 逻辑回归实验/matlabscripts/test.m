       clear all;%%��ʾ������б���������ֵ%%

       t=0:pi/360:2*pi;

       x=sin(t);

       y=cos(t);

       z=2*x.^2+y.^2;

       scatter3(1,0,1);
       
       hold on
     
       plot3(x,y,z,'Color','r','LineWidth',2);

      %%��ά����������ͱ��������%%

      xlabel('x');

      ylabel('y');

      zlabel('z');

      title('��ά����ͼ');
      
      view(40,35)
