function fill_plot(x,y,stdx,stdy,c)

if size(x,1) > 1 && size(x,2) == 1                % If not given in row arrays
    x = x';
    y = y';    
    stdx = stdx';
    stdy = stdy';
end
fill([x+stdx, fliplr(x-stdx)],[y+stdy, fliplr(y-stdy)],c,'FaceAlpha',0.2,'EdgeColor','none');
