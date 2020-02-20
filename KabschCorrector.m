% Correct tracks to average position

function [tracks_corrected,a_all,d_avg,N] = KabschCorrector(tracks,number_Frames)
%matlabpool(6)
%tic
d_avg = zeros(number_Frames,1);
N = zeros(number_Frames,1);
a_all = zeros(number_Frames,5);
tracks_corrected = tracks;
ntracks = max(tracks(:,end));
tracks_avg = zeros(ntracks,3);
for id = 1:ntracks
    tracks_avg(id,1:2) = mean(tracks(tracks(:,end) == id,1:2),1);
    tracks_avg(id,3) = id;
end

% frame_start_i = min(tracks(:,end-1));
% for frame_start = frame_start_i:number_Frames
%     if dim == 3
%         a_all(frame_start,:) = [0 0 0 0 0 0 1 0 0 0];
%     else
%         a_all(frame_start,:) = [0 0 0 0 0 1 0 0];
%     end
%     tracks_old = tracks(tracks(:,end-1)==frame_start,:);
%     frame_end = frame_start + dt;
%     if frame_end > number_Frames
%         frame_end = number_Frames;
%     end
    for i = 1:number_Frames
        %N_part = 60; % maximum number of particles to consider for approximating affine transform (was 10 for initial correction)
        tracks_i = tracks(tracks(:,end-1)==i,:);
        [~,ind_avg,ind_i] = intersect(tracks_avg(:,end),tracks_i(:,end));
        avg_pos = tracks_avg(ind_avg,1:2);
        i_pos = tracks_i(ind_i,1:2);
        intens = tracks_i(ind_i,3);
        [~,I] = sort(intens);
        I = flipud(I);
        N(i) = length(I);
        [U,r,d_avg(i)] = Kabsch(i_pos(I,:)',avg_pos(I,:)',intens');
        
        %[U,r,d_avg(i)] = Kabsch(i_pos(I,:)',avg_pos(I,:)'); This line can
        %be used instead if you don't want to weight by intensity
        
        %else
        %    [U,r,d_avg(i)] = Kabsch(i_pos(I(1:N_part),:)',avg_pos(I(1:N_part),:)');
        %end
        %e_out = EulerParams(U);
        %theta_out = 2*acos(e_out(1));
        %axis_out = e_out(2:4)/sin(theta_out/2);
        theta_out = asin(U(2,1));
        a_all(i,1:3) = [r(1:2)' theta_out];
        i_pos_2correct = tracks(tracks(:,end-1)==i,1:2);
        num = size(i_pos_2correct,1);
        i_pos_shift = i_pos_2correct + repmat(r',num,1);
        i_pos_shift_c = (i_pos_shift'*ones(num,1)/num)';
        a_all(i,4:5) = i_pos_shift_c'; % Centroid of translated(!) image about which coordinates are rotated
        i_pos_axis = i_pos_shift - repmat(i_pos_shift_c,num,1);
        i_pos_axis(:,2) = -i_pos_axis(:,2);
        i_pos_rot = zeros(num,2);
        for p = 1:num
            r_new = i_pos_axis(p,:)';
            i_pos_rot(p,:) = (U*r_new)';
        end
        i_pos_rot(:,2) = -i_pos_rot(:,2);
        i_pos_final = i_pos_rot + repmat(i_pos_shift_c,num,1);
        tracks_corrected(tracks_corrected(:,end-1)==i,1:2) = i_pos_final;
    end
%end
clear tracks_old tracks_new intens I a_out

%matlabpool close