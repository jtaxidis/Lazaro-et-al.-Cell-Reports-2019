function plot_hist_signif(p,x1,x2,y1,y2,incrx,incry)

mx = mean([x1 x2]);
my = mean([y1 y2]);

if p <= 0.05
    line([x1,x2]+incrx , [y1 y2]+incry,'Color','k','Linewidth',1.5);
    plot(mx + 2*incrx , my+2*incry,'*m','MarkerSize',10);
end

