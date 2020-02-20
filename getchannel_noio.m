function out = getchannel_noio(pthname,fname,savefolder,dual_ch)

% This is all a segment from Dan Winkowski's code

% disp('Select .raw file to segment into images');
% 
% % SELECT RAW IMAGE TO SEGMENT - BIG Thnx to PW for advice
% if ismac == 0
%     [fname, pthname] = uigetfile([rawfolder '*.raw'], 'Select a RAW DataFile?');
%     fullpathIMG = [pthname fname];
% else
%     [fname, pthname] = uigetfile('/Users/zbowen/Documents/research/VolScan0001/TSeries0007/*.raw','Select a RAW DataFile?');
%     fullpathIMG = [pthname fname];
% end

fullpathIMG = [pthname fname];

% Create directories within savefolder
if dual_ch
    if ~ismac
        mkdir([savefolder '\redChan\']);
        mkdir([savefolder '\redChan\Analysis\']);
        mkdir([savefolder '\greenChan\']);
        mkdir([savefolder '\greenChan\Analysis\']);
    else
        mkdir([savefolder '/redChan/']);
        mkdir([savefolder '/redChan/Analysis/']);
        mkdir([savefolder '/greenChan/']);
        mkdir([savefolder '/greenChan/Analysis/']);
    end
end

% Read associated timing.txt file for timing index to image sequences
TxtFname = [pthname 'timing.txt'];
fileID = fopen(TxtFname);
formatSpec = '%f';
TimingVals = fscanf(fileID,formatSpec);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For this testing, got rid of a section of code right here %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read associated XML file for image dimensions
fullpathXML = [pthname 'Experiment.xml'];
tt = ParseXML(fullpathXML);
ChanAEnabled = str2double(tt.Children(30).Attributes(1).Value);
ChanBEnabled = str2double(tt.Children(30).Attributes(2).Value);
ChanCEnabled = str2double(tt.Children(30).Attributes(3).Value);
ChanDEnabled = str2double(tt.Children(30).Attributes(4).Value);
numChannels = ChanAEnabled + ChanBEnabled +ChanCEnabled +ChanDEnabled;
numTimePts = str2double(tt.Children(36).Attributes(5).Value);
numImages = numChannels * numTimePts;
dimXmicrons = str2double(tt.Children(22).Children(2).Children(2).Attributes(2).Value) *1000;
dimYmicrons = str2double(tt.Children(22).Children(2).Children(2).Attributes(3).Value) * 1000;
dimXpixels = str2double(tt.Children(28).Attributes(16).Value);
dimYpixels = str2double(tt.Children(28).Attributes(17).Value);
micPerPixelX = dimXmicrons / dimXpixels;
micPerPixelY = dimYmicrons / dimYpixels;
% approxNeuronDiamRangeUM = [5 14]; % assumption of neuron diameter to be between 6 & 14 microns
% approxNeuronHWRangePx = round((approxNeuronDiamRangeUM ./ micPerPixelX)/2); %divide by 2 becasue its a radius
%str = input('Select one or more channel (example:2 3) :','s');
% if dual_ch == 4
%     str = '2 3';
% else
%     str = '1';
% end
%selectedChan = str2num(str);
% numSelected = length(selectedChan);
numTotImages = str2double( tt.Children(36).Attributes(5).Value );
numTotTZImages = numTotImages * numChannels;

numBytes = 2; % assumes 16-bit
segmentSize = (dimXpixels * dimYpixels * numBytes); %in case multiple channels are enabled (ChanB is green)
% skipIM = (numChannels-numSelected) * segmentSize;  %in case multiple channels are enabled (ChanB is green)
out.field_width = numel(num2str(length(TimingVals)))+1;
clear tt

% Read RAW file into workspace
if numChannels == 1
    fh = fopen(fullpathIMG); % or whatever the filename is
    status = fseek(fh,0,'bof'); %offset to Green Chan (ChanB)
    pos = ftell(fh);
    disp('pos : ',num2str(pos))
    IMG = reshape(fread(fh,inf,'uint16=>uint16'),[dimXpixels, dimYpixels, numTimePts]);
    IMG = permute(IMG, [2 1 3]);
    fclose(fh);
    mnIMG = squeeze(mean(IMG(:,:,:),3));  % bc the 1st iteration in ThorImage is "split screen"
    telapsedReading = toc(tstart);
    disp(['Time Elapsed for Reading was: ', num2str(telapsedReading/60), ' minutes'])
else
    % preallocations
    %a = zeros(4,1);
    %a(selectedChan) = 1; %channel selection vector
    if dual_ch
        %file(length(TimingVals)).channelB = [];
        %out.IMG_B = zeros(dimXpixels,dimYpixels,length(TimingVals));
        %file(length(TimingVals)).channelC = [];
        out.IMG_B = zeros(dimXpixels,dimYpixels,length(TimingVals));
        
        out.IMG_D = zeros(dimXpixels,dimYpixels,length(TimingVals));
    else
        %file(length(TimingVals)).channelA = [];
        out.IMG_A = zeros(dimXpixels,dimYpixels,length(TimingVals));
    end
    
       %% read entire raw file first
    
%    fh = fopen(fullpathIMG); % or whatever the filename is
%    status = fseek(fh,0,'bof') %offset to Green Chan (ChanB)
%     pos = ftell(fh);
%     disp('pos : ',num2str(pos))
    fh = fopen(fullpathIMG); % opens the file
    IMGIMG=fread(fh,inf,'uint16=>uint16');
    IMGread = reshape(IMGIMG, [dimXpixels, dimYpixels,4*length(TimingVals)]);
    clear IMGIMG
    
    IMGmultiCh = IMGread;
    clear IMGread
    fclose(fh);
    for i=1:length(TimingVals)
        kB= 4*(i-1)+2;
        kD=4*(i-1)+4; % +4 for IMG_D
        out.IMG_B(:,:,i)=IMGmultiCh(:,:,kB);
        
        out.IMG_D(:,:,i)=IMGmultiCh(:,:,kD);
    end
    out.IMG_B(:,:,:) = permute(out.IMG_B(:,:,:), [2 1 3]);
    out.IMG_D(:,:,:) = permute(out.IMG_D(:,:,:), [2 1 3]);
    out.IMG_C=out.IMG_D;
%     fh = fopen(fullpathIMG); % opens the file
% 
% % loops through time, only reads selected channels
% % due to preallocation, there will still be an IMG for each channel, but
% % the unselected channels will remain empty matrices
%     disp(['Parsing ' num2str(length(TimingVals)) ' images...']);
%     if dual_ch
%         for i = 1:length(TimingVals)
%         %disp(['Image ' num2str(i) '...']);
%             fseek(fh,segmentSize,'cof');
%             %file(i).channelB = fread(fh , segmentSize/2 , 'uint16' , 'l'); 
%             %fseek(fh,segmentSize,'cof');
%             
%             out.IMG_B(:,:,i) = reshape(fread(fh , segmentSize/2 , 'uint16' , 'l'),[dimXpixels, dimYpixels]);
%             out.IMG_B(:,:,i) = permute(out.IMG_B(:,:,i), [2 1 3]);
%             %file(i).channelC = fread(fh , segmentSize/2 , 'uint16' , 'l');
%             out.IMG_C(:,:,i) = reshape(fread(fh , segmentSize/2 , 'uint16' , 'l'),[dimXpixels, dimYpixels]);
%             out.IMG_C(:,:,i) = permute(out.IMG_C(:,:,i), [2 1 3]);
%             fseek(fh,segmentSize,'cof');
%             %clear file
%         end
%     else
%         for i = 1:length(TimingVals)
%             out.IMG_A(:,:,i) = reshape(fread(fh , segmentSize/2 , 'uint16' , 'l'),[dimXpixels, dimYpixels]);
%             out.IMG_A(:,:,i) = permute(out.IMG_B(:,:,i), [2 1 3]);
%             fseek(fh,3*segmentSize,'cof');
%         end
%     end
% end
% 
% 
% % % this loop plays a movie going through all of the frames
% % figure
% % for i=1:size(IMG_B,3)
% % imagesc(IMG_B(:,:,i),[0 max(IMG_B(:))])
% % axis off
% % drawnow
% % pause(1/60) %denominator is the FPS
 end