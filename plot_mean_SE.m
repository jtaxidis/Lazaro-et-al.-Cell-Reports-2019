function plot_mean_SE(x,data,pvalue,cols)

if isempty(x)
    x = 1:length(data);
end

if isscalar(cols)
    col = cols;
    for i = 1:length(data)
        cols(i) = col;
    end
end

hold on;
for i = 1:length(data)
    M = nanmean(data{i});
    S = nanstd(data{i})/sqrt(length(data{i}));
    
    bar(x(i),M,0.3,'EdgeColor',cols(i),'FaceColor','none','LineWidth',2);
%     plot(i+rand(length(data{i}),1)*0.2-0.1,data{i},'o','Color',cols(i));
    plot(x(i)*ones(length(data{i}),1),data{i},'o','Color',cols(i),'MarkerSize',4);
    errorbar(x(i),M,S,'k','Linewidth',2.5);
end

if ~isempty(pvalue)
    for i = 1:length(data)-1
        for j = i+1 : length(data)
            plot_signif(pvalue(i,j),x(i),x(j),max([mean(data{i}) mean(data{j})]),mean(data{i})*0.05);
        end
    end
end


function mb = plot_signif(p,loc1,loc2,mbold,incr)

if p <= 0.05
    mb = mbold + incr*4;
    ml = mean([loc1,loc2]);
    
    line([loc1,loc1] , [mb, mb+incr],'Color','k');
    line([loc2,loc2] , [mb, mb+incr],'Color','k');
    line([loc1,loc2] , [mb+incr, mb+incr],'Color','k');
    
    if p <= 0.05 && p > 0.01
        plot(ml , mb+incr*2,'*k','MarkerSize',5);
    elseif p <= 0.01 && p > 0.001
        plot(ml + [-ml/50,ml/50] , [mb+incr*2, mb+incr*2] , '*k','MarkerSize',5);
    elseif p <= 0.001
        plot(ml + [-ml/25,0,ml/25] , [mb+incr*2, mb+incr*2, mb+incr*2] , '*k','MarkerSize',5);
    end
else
    mb = mbold;
end

