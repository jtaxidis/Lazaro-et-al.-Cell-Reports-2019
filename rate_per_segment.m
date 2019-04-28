function R = rate_per_segment(S,t1,t2,binlen)

btime = t1 : binlen : t2;                                                   % Make bin times within the segment
lb = length(btime);
if lb > 0
    btime(end) =  t2;                                                           % (In case the last bin was a bit shorter)
end
R = histcounts(S,btime);
R = R/binlen;