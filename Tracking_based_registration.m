%% Master code for tracking based registration

gclc
clear
%% Selecting Data 

dual_ch = input('How many channels is your data? 1 or 4? ');
red_tr = input('Would you like to track the red or green channel? 1 for red, 0 for green. '); % Set to 1 if the raw file includes both channels; 0 if they are separate raw files
disp('Select .raw file to segment into images');
if ~ismac
    [fname, pthname] = uigetfile([pwd '\*.raw'], 'Select a RAW DataFile');
    fullpathIMG = [pthname fname];
else
    [fname, pthname] = uigetfile([pwd '/*.raw'],'Select a RAW DataFile');
    fullpathIMG = [pthname fname];
end
disp('Select/create folder to save analysis');
savefolder = uigetdir(pwd, 'Select an Analysis Folder');
%% parsing images
% This section reads the raw file from the location selected in the previous section and
% parses each channel into frames. 

%parpool(4);
tic
out = getchannel_noio(pthname,fname,savefolder,dual_ch); % reading the XML file that 
% stores two-photon imaging metadata 
field_width=out.field_width;

if dual_ch == 1
    IMG_tr = out.IMG_A;
else
    if red_tr
        IMG_tr = out.IMG_C;
        out.IMG_C = [];
        IMG_sig = out.IMG_B;
        out.IMG_B = [];
    else
        IMG_tr = out.IMG_B;
        out.IMG_B = [];
        IMG_sig = out.IMG_C;
        out.IMG_C = [];
    end
end
clear out
toc
[Dx,Dy,Dz]= size(IMG_tr);
%%
%refine = input('Is this a refinement of registered images? 1 for yes, 0 for no. ');
if red_tr == 1
%    if refine
%        file_short='red_reg';
%        if ~ismac
%            pathname = [savefolder '\redChan\Analysis\Registration\'];
%        else
%            pathname = [savefolder '/redChan/Analysis/Registration/'];
%        end
%    else
%        file_short='red';
        if ~ismac
            pathname = [savefolder '\redChan\'];
        else
            pathname = [savefolder '/redChan/'];
        end
%    end
else
%    if refine
%        file_short='green_reg';
%        if ~ismac
%            pathname = [savefolder '\greenChan\Analysis\Registration\'];
%        else
%            pathname = [savefolder '/greenChan/Analysis/Registration/'];
%        end
%    else
%        file_short='green';
        if ~ismac
            pathname = [savefolder '\greenChan\'];
        else
            pathname = [savefolder '/greenChan/'];
        end
%    end
end  % Image file name, not including numbers


frame_start = input('What is the first frame you would like to track/analyze? ');
frame_end = input('What is the last frame you would like to track/analyze? ');
number_Frames = frame_end-frame_start+1;
disp('Would you like to generate and save position overlay images? A test image of the first frame will be saved regardless.')
plot_on = input('1 for yes, 0 for no. ');%0
if plot_on == 1
    frame_start_plot = 1;%30input('What is the first frame you would like to save? ');
    plot_interval = 5;%input('What is the interval between saved extraction images (e.g., 1 will save every image)? ');
end
if ~ismac
    mkdir([pathname 'Analysis\Extraction\']);
    savefolder_extr = [pathname 'Analysis\Extraction\'];
else
    mkdir([pathname 'Analysis/Extraction/']);
    savefolder_extr = [pathname 'Analysis/Extraction/'];
end

%% setting initial settings for tracking

loop = 1;
while loop
    lower_bpass = input('What is the typical size (in pixels) of noise in your images? '); % 4 for Dan's data
    higher_bpass = input('What is the size (in pixels) of features (i.e., cells) you would like to find? '); % Should be odd number; 
    N_user = input('How many cells would you like to extract? ');
    %threshold = input('What intensity threshold would you like to impose (this is used for identifying significant maxima)? '); %275;
    [pos_lst_test,N_p] = CellExtractor_noio(IMG_tr(:,:,frame_start),frame_start,lower_bpass,higher_bpass,N_user,1,savefolder_extr,1);
    disp('The bandpassed image and extracted features (cells) are shown for the first frame.');
    loop = input('Would you like to modify any settings? 1 for yes, 0 for no. ');
    close all
end
%% Cell Extract
pos_lst_cell = cell(number_Frames,1);
%N = input('How many cores would you like to use? ');
%parpool(N);
tic
for i = frame_start:frame_end %parfor changed
    pos_lst_cell{i} = CellExtractor_noio(IMG_tr(:,:,i),i,lower_bpass,higher_bpass,N_user,(plot_on) && (mod(i-frame_start_plot,plot_interval)==0),savefolder_extr,0);
close all
end
toc
%delete(gcp);
pos_lst = cell2mat(pos_lst_cell);
clear pos_lst_cell
file_short='red_reg';
%%
param.quiet = 1;
param.dim = 2;
plot_on_tr = input('Would you like to generate and save track images? 1 for yes, 0 for no. ');%0
if plot_on_tr
    if plot_on
        keep = input('Would you like to keep the same settings for frequency of saving extraction images? 1 for yes, 0 for no. ');
        if keep
            frame_start_plot_tr = frame_start_plot;
            plot_interval_tr = plot_interval;
        end
    else
        frame_start_plot_tr = input('What is the first frame you would like to save? ');
        plot_interval_tr = input('What is the interval between saved extraction images (e.g., 1 will save every image)? ');    
    end
end
check_tr = 0;
loop_tr = 1;
while loop_tr || check_tr
    max_disp = input('What is the furthest distance (in pixels) that features can move between frames? ' ); %10;
    param.good = input('How long (in frames) must a track last to be valid? ' ); %2500;
    param.mem = input('How many frames is a feature allowed to be missed without starting a new track? ' ); %max(pos_lst(:,end)); % I have memory maxed out
    tic
    tracks = track(pos_lst,max_disp,param);
    toc
    ntracks = max(tracks(:,end));
    disp([num2str(ntracks) ' total tracks were found.']);
    Nt = hist(tracks(:,end-1),1:number_Frames);
    figure(1)
    plot(1:number_Frames,Nt);
    xlabel('Frame Number');
    ylabel('No. of Active Tracks');
    check_tr = input('Would you like to try new track settings before continuing? 1 for yes, 0 for no. ');
    if check_tr
        continue;
    end

    if plot_on_tr
        if ~ismac
            savefolder_tr = [pathname 'Analysis\Tracking\'];
        else
            savefolder_tr = [pathname 'Analysis/Tracking/'];
        end
        col = repmat(jet(8),ceil(ntracks/8),1);
        particle_size = higher_bpass;
        [circle_x,circle_y] = jcircle(particle_size/2);
        for t = frame_start_plot_tr:plot_interval_tr:number_Frames
            filename=strcat(pathname,sprintf(['%s%0' num2str(field_width) 'u.tif'],file_short,t));
            a = double(imread(filename));
            h = figure(3);
            figure(h);
            colormap gray
            imagesc(a);
            truesize;
            hold on
            tracks_t = tracks(tracks(:,end-1) == t,:);
            for j = 1:size(tracks_t,1)
                plot(circle_x+tracks_t(j,1),circle_y+tracks_t(j,2),'Color',col(tracks_t(j,end),:),'LineWidth',1);
            end
            axis off
            axis equal
            set(gca,'units','pixels') % set the axes units to pixels
            x = get(gca,'position'); % get the position of the axes
            set(gcf,'units','pixels'); % set the figure units to pixels
            y = get(gcf,'position'); % get the figure position
            set(gcf,'position',[y(1) y(2) x(3) x(4)]);% set the position of the figure to the length and width of the axes
            set(gca,'units','normalized','position',[0 0 1 1]); % set the axes units to pixels
            set(gcf, 'PaperPositionMode', 'auto');
            savename=strcat(savefolder_tr,sprintf(['tracks_%0' num2str(field_width) 'u'],t),'.tif');
            print(gcf,savename,'-dtiff');
            hold off;
        end
        loop_tr = input('From observing the figures, would you like to adjust any track settings? 1 for yes, 0 for no. ');
    else
        loop_tr = 0;
    end
end
ntracks = max(tracks(:,end));
tracks_avg = zeros(ntracks,param.dim+1);
for id = 1:ntracks
    tracks_avg(id,1:param.dim) = mean(tracks(tracks(:,end) == id,1:param.dim),1);
    tracks_avg(id,end) = id;
end
%%
% Generate corrected tracks

d_avg_raw = rms_disp(tracks,tracks_avg,frame_start,frame_end,param.dim);
tic
[tracks_corrected,a_all,d_avg,Ntr_corr] = KabschCorrector(tracks,number_Frames);
toc
%%

% !! This could be rewritten as a function

% Apply transforms to images
% First translate, then rotate about the centroid of common cells
%if refine
%    mkdir([savefolder '\redChan\Analysis\Registration\Analysis\Registration\']);
%    mkdir([savefolder '\greenChan\Analysis\Registration\Analysis\Registration\']);
%else
if ~ismac
    mkdir([savefolder '\greenChan\Analysis\Registration\']);
    mkdir([savefolder '\redChan\Analysis\Registration\']);
else
    mkdir([savefolder '/greenChan/Analysis/Registration/']);
    mkdir([savefolder '/redChan/Analysis/Registration/']);
end
%end
%red_avg_reg = 0;
%green_avg_reg = 0;
IMG_tr_reg = zeros(size(IMG_tr,1),size(IMG_tr,2),number_Frames);
IMG_sig_reg = IMG_tr_reg;
pad_d = max(round(max(abs(a_all(:,1:2))) + size(IMG_tr,1)*max(a_all(:,3))));
%parpool(N);
tic
for i = 1:number_Frames %parfor changed
    
    
    IMG_tr_reg(:,:,i) = correct_image(IMG_tr(:,:,i),a_all(i,:),pad_d);
    IMG_sig_reg(:,:,i) = correct_image(IMG_sig(:,:,i),a_all(i,:),pad_d);
    %if refine
    %    im_red = imread([savefolder '\redChan\Analysis\Registration\red_reg' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
    %    im_green = imread([savefolder '\greenChan\Analysis\Registration\green_reg' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
    %else
        
    
    %red_avg_reg = red_avg_reg + double(im_red_reg);
    %green_avg_reg = green_avg_reg + double(im_green_reg);
    %if refine
    %    imwrite(IMG_A_reg,[savefolder '\redChan\Analysis\Registration\Analysis\Registration\red_reg' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
    %    imwrite(IMG_B_reg,[savefolder '\greenChan\Analysis\Registration\Analysis\Registration\green_reg' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
    %else
    %    imwrite(IMG_A_reg,[savefolder '\redChan\Analysis\Registration\red_reg' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
    %    imwrite(IMG_B_reg,[savefolder '\greenChan\Analysis\Registration\green_reg' num2str(i,['%0' num2str(field_width) 'd']) '.tif']);
    %end
end
IMG_tr_reg_avg = mean(double(IMG_tr_reg),3);
IMG_sig_reg_avg = mean(double(IMG_sig_reg),3);
% red_avg_reg = red_avg_reg/number_Frames; 
% green_avg_reg = green_avg_reg/number_Frames; 
toc
%delete(gcp);
%%
% Extraction of cells from the average registered red image
donut_check = 1;
close all
while donut_check
    lower_bpass_avg = input('What is the typical size (in pixels) of noise in your images? '); % 4 for Dan's data
    higher_bpass_avg = input('What is the size (in pixels) of features (i.e., cells) you would like to find? '); % Should be odd number; 
    b_avg=bpass(IMG_C_reg_avg,lower_bpass_avg,higher_bpass_avg);
    figure(1);
    colormap gray
    %imagesc(red_avg_reg);
    imagesc(IMG_C_reg_avg);
    truesize;
    set(gca,'units','pixels') % set the axes units to pixels
    x = get(gca,'position'); % get the position of the axes
    set(gcf,'units','pixels'); % set the figure units to pixels
    y = get(gcf,'position'); % get the figure position
    set(gcf,'position',[y(1) y(2) x(3) x(4)]);% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0 0 1 1]); % set the axes units to pixels
    set(gcf, 'PaperPositionMode', 'auto');
    figure(2);
    imagesc(b_avg)
    truesize;
    set(gca,'units','pixels') % set the axes units to pixels
    x = get(gca,'position'); % get the position of the axes
    set(gcf,'units','pixels'); % set the figure units to pixels
    y = get(gcf,'position'); % get the figure position
    set(gcf,'position',[y(1) y(2) x(3) x(4)]);% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0 0 1 1]); % set the axes units to pixels
    set(gcf, 'PaperPositionMode', 'auto');
    thresh_avg = input('What threshold value would you like to impose? ');
    pk_avg = pkfnd(b_avg,thresh_avg,higher_bpass_avg);
    cnt_avg = cntrd(IMG_C_reg_avg,b_avg,pk_avg,higher_bpass_avg);
    if (~isempty(cnt_avg))
        [circle_x,circle_y] = jcircle(higher_bpass_avg/2);
        figure(1);
        hold on
        for j = 1:size(cnt_avg,1)
            plot(circle_x+cnt_avg(j,1),circle_y+cnt_avg(j,2),'b-','LineWidth',1);
        end
        hold off
        pos_lst_avg = zeros(length(cnt_avg(:,1)),5);
        pos_lst_avg(:,1:5)=cnt_avg;
        donut_check = input('Would you like to try new input parameters? 1 for yes, 0 for no. ');
    else
        disp('No cells found!');
    end
end
clear cnt_avg