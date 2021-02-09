function [Avisoft_times,Open_ephys_times] = avisoftSync(params_file)
%% load the files and create the variables
run(params_file);
data = load_open_ephys_binary([exp_path,'/experiment1/recording1/structure.oebin'],'continuous',1,'mmap');
channel_number = data.Header.num_channels;
aud_op = data.Data.Data.mapped(channel_number-(8-mic),:);
[avi_file,dir] = uigetfile;
[aud_avi] = audioread([dir,avi_file]);

%% Find start and end point in both variables
avi_first_section = aud_avi(1:300000);
avi_start = find(abs(avi_first_section) > iqr(avi_first_section) * 10,1);
avi_last_section = aud_avi(end-300000:end);
avi_end = find(abs(avi_last_section) > iqr(avi_last_section) * 10);
avi_end = length(aud_avi) - 300000 + avi_end(1);

open_first_section = aud_op(1:45000);
open_start = find(abs(open_first_section) > iqr(open_first_section) * 20,1);
open_last_section = aud_op(end-45000:end);
open_end = find(abs(open_last_section) > iqr(open_last_section) * 20,1);
open_end = length(aud_op) - 45000 + open_end;

length_open = open_end - open_start;
time_total = length_open/SF;
true_avi_SF = (avi_end - avi_start)/time_total;

time_vector_open = [1:(open_end - open_start)]/SF;
time_vector_avi = [1:(avi_end - avi_start)]/true_avi_SF;







%% create output structures
Avisoft_times.start = avi_start;
Avisoft_times.end = avi_end + 1;
Avisoft_times.timestamps = time_vector_avi;
Avisoft_times.SF = true_avi_SF;
Open_ephys_times.start = open_start + 1;
Open_ephys_times.end = open_end;
Open_ephys_times.timestamps = time_vector_open;
end

