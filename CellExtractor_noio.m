function [pos_lst,N_p] = CellExtractor_noio(IMG_tr,frame,lower_bpass,higher_bpass,N_user,plot_on,savefolder,test_im)

particle_size=higher_bpass;

%filename=strcat(pathname,sprintf(['%s%0' num2str(field_width) 'u.tif'],file_short,frame));
a = double(IMG_tr);
%a1 = padarray(a,[higher_bpass,higher_bpass],'both');
    
b=bpass(a,lower_bpass,higher_bpass);
if test_im
    figure(1);
    imagesc(b);
    truesize;
    set(gca,'units','pixels') % set the axes units to pixels
    x = get(gca,'position'); % get the position of the axes
    set(gcf,'units','pixels'); % set the figure units to pixels
    y = get(gcf,'position'); % get the figure position
    set(gcf,'position',[y(1) y(2) x(3) x(4)]);% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0 0 1 1]); % set the axes units to pixels
    set(gcf, 'PaperPositionMode', 'auto');
end
p = 99.99;
pk = [];
step = -0.01;
diff = -1;
N_p = 0;

while (size(pk,1) ~= N_user)
    N_p = N_p + 1;
    threshold = prctile(b(:),p);
    pk = pkfnd(b,threshold,particle_size);
    diff_new = sign(size(pk,1)-N_user);
    if diff_new ~= diff;
        step = -0.1*step;
        diff = diff_new;
    end
    p = p + step;
end

if (~isempty(pk)) 
    cnt = cntrd(a,b,pk,particle_size);
    if (~isempty(cnt))
        if plot_on
            [circle_x,circle_y] = jcircle(particle_size/2);
            h = figure(2);
            figure(h);
            colormap gray
            imagesc(a);
            truesize;
            hold on
            for j = 1:size(cnt,1)
                plot(circle_x+cnt(j,1),circle_y+cnt(j,2),'b-','LineWidth',1);
            end
        end
        pos_lst = zeros(length(cnt(:,1)),6);
        pos_lst(:,1:5)=cnt(:,1:5);
        pos_lst(:,6)=frame;   
    end
end
if plot_on
    axis off
    axis equal
    set(gca,'units','pixels') % set the axes units to pixels
    x = get(gca,'position'); % get the position of the axes
    set(gcf,'units','pixels'); % set the figure units to pixels
    y = get(gcf,'position'); % get the figure position
    set(gcf,'position',[y(1) y(2) x(3) x(4)]);% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0 0 1 1]); % set the axes units to pixels
    set(gcf, 'PaperPositionMode', 'auto');
    savename=strcat(savefolder,sprintf('overlay_%05u',frame),'.tif');
    print(gcf,savename,'-dtiff');
    hold off;
end