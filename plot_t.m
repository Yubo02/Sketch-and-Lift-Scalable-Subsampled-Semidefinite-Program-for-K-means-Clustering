


x=    tim1_p(2,:);
y1 =  tim1_p(1,:);
y2 =  tim2_p(1,:);
y3 =  tim3_p(1,:);
y4 =  tim4_p(1,:);
y5=   tim5_p(1,:);



figure;
set(gcf,'unit','centimeters','position',[10 6 20 15]);
set(gcf,'Color',[0.9 0.9 0.9]);

axis([min(x) max(x) 0 10]);

grid on;

figure_FontSize=8;
xlabel({'p'},'FontSize',figure_FontSize);
ylabel({'Time cost (s)', ' ', ' '});

set(get(gca,'XLabel'),'FontSize',12,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',12,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(gca,'Position',[.13 .17 .80 .74]);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
title('Time cost when p changes');

hold on
p=plot(x,y1,'-*',x,y2,'--^',x,y3,':d',x,y4,'-<',x,y5,'-s');
 p(1).LineWidth = 2;
 p(1).Color='#A2142F';
 p(2).LineWidth = 2;
 p(2).Color='#7E2F8E';
 p(3).LineWidth = 2;
 p(3).Color='#000000';
 p(4).LineWidth = 2;
 p(4).Color='#0072BD';
 p(5).Color='m';
 p(5).LineWidth = 2;
hold off


hold off


legend(p,'M_1','M_2','M_3','M_4','M_5','Location','northeast','Orientation','vertical');




