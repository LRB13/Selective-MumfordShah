function [gd,mask,cols,rows] = geodist(Im,seg,cols,rows)

Im = double(Im);
Im = Im-min(Im(:)); Im=Im./max(Im(:));

[n,m] = size(Im);

if nargin < 2 %If only image is input, get user input
    figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    [mask,cols,rows] = roipoly(Im);
    mask = double(mask);
    if norm(mask,2) ==0
        for i=1:size(cols,1)
        mask(round(rows(i)),round(cols(i)))=1;
        end
    end
else %if n=2,3 or 4
    seg(seg>=.5) = 1;
    seg(seg<=0.5) = 0;
    mask = seg;
end

if nargin ==2
    cols=-1;
    rows=-1;
end

if nargin > 2
        mask = seg+roipoly(Im,cols,rows);
    mask(mask>.5) = 1; mask(mask<=.5) = 0;
end


%make mask bigger for very slightly quicker performance
if sum(mask(:))==1
    [x,y] = find(mask==1);
    mask(x+1,y) = 1; mask(x-1,y) = 1;
    mask(x,y+1) = 1; mask(x,y-1) = 1;
    mask(x+1,y+1) = 1;
end

%%%%%%%%%%%%%%%%%% optional smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ims = imgaussfilt(Im,2);
ims = Im;
%ims = anisodenoise(Im);
ims = imtgvsmooth(Im,.05,.025,10);
[gx,gy] = gradient(ims);
grad = gx.^2 + gy.^2;
%%%

eucdist = timesweep(ones(n,m),mask);

beta = 1000;
epsi = 1e-3;
theta=0; %in paper original paper, theta = 0.1.

f = epsi + 1000.*grad + theta.*eucdist;
tic
D = timesweep(f,mask);
toc
D = D - min(D(:));
gd=D./max(D(:));
end