function d_avg = rms_disp(tracks, tracks_avg, frame_start, frame_end, dim)

d_avg = zeros(frame_end-frame_start+1,1);
for i = frame_start:frame_end
    new_pos = tracks(tracks(:,end-1)==i,:);
    ids = intersect(new_pos(:,end),tracks_avg(:,end));
    avg_ids = tracks_avg(ismember(tracks_avg(:,end),ids),1:dim);
    new_ids = new_pos(ismember(new_pos(:,end),ids),1:dim);
    d_avg(i) = sqrt(mean(sum((new_ids - avg_ids).^2,2)));
end