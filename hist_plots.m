clear all;
a_data
% Histogram plots ---------------------------------------------------------
% a_data
nbins=9;
hist(a_data,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');%
title("Histogram Values for 'a' Data",'fontname', 'Times', 'fontsize',20);
xlabel('','fontname', 'Times', 'fontsize',20);
ylabel('','fontname', 'Times', 'fontsize',20);
grid on;
set(gca,'FontSize',16);
print -dpng histogram_Values_a_10000.png
pause;
clf;
