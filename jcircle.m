function [outx,outy]=jcircle(r)
cnt=0;
for theta=0:.5:2*pi+.5
    cnt=cnt+1;
    outx(cnt)=r*cos(theta);
    outy(cnt)=r*sin(theta); 
end

