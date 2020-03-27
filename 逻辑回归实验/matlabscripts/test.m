       clear all;%%表示清空所有变量及变量值%%

       t=0:pi/360:2*pi;

       x=sin(t);

       y=cos(t);

       z=2*x.^2+y.^2;

       scatter3(1,0,1);
       
       hold on
     
       plot3(x,y,z,'Color','r','LineWidth',2);

      %%三维曲线坐标轴和标题的设置%%

      xlabel('x');

      ylabel('y');

      zlabel('z');

      title('三维曲线图');
      
      view(40,35)
