function [p,testtype] = significance(dist1,dist2,type)

if strcmp(type,'unequal')
    ttype = 'both';
elseif strcmp(type,'larger')
    ttype = 'right';
elseif strcmp(type,'smaller')
    ttype = 'left';
end

p = 1;
testtype = 0;


h1 = 1;
h2 = 1;
if sum(~isnan(dist1)) >= 4  
    h1 = lillietest(dist1,0.05,'norm');
end

if length(dist2) == 1
    if h1 == 1
        p = signtest(dist1,dist2,'alpha',0.05,'tail',ttype);
        testtype = 1;
    else
        [~,p] = ttest(dist1,dist2,0.05,ttype);
        testtype = 2;
    end
    
else
    if sum(~isnan(dist2)) >= 4
        h2 = lillietest(dist2,0.05,'norm');
    end
    
    if length(dist1) > 1 && length(dist2) > 1
        if any([h1 h2] == 1)
            p = ranksum(dist1,dist2,'alpha',0.05,'tail',ttype);
            testtype = 3;
        else
            [~,p] = ttest2(dist1,dist2,[],ttype,'unequal');
            testtype = 4;
        end
    end
end
