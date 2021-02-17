function [Avisoft_times] = avisoftSync(params_file)
%% load the files and create the variables
run(params_file);
data = load_open_ephys_binary([exp_path,'/experiment1/recording1/structure.oebin'],'continuous',1,'mmap');
channel_number = data.Header.num_channels;
aud_op = data.Data.Data.mapped(channel_number-(8-mic),:);
cd(exp_path);
files = dir('*.wav');
[aud_avi,avi_soft_sf] = audioread([files.folder,'/',files.name]);

%% Find start and end point in both variables
avi_end_frame = 61 * avi_soft_sf;
avi_first_section = aud_avi(1:300000);
avi_start = find(abs(avi_first_section) > iqr(avi_first_section) * 10,1);
avi_last_section = aud_avi(end-avi_end_frame:end);
avi_end = find(abs(avi_last_section) > iqr(avi_last_section) * 10);
avi_end = length(aud_avi) - avi_end_frame + avi_end(1);


open_end_frame = 61*30000;
open_first_section = aud_op(1:45000);
open_start = find(abs(open_first_section) > iqr(open_first_section) * 20,1);
open_last_section = aud_op(end-open_end_frame:end);
open_end = find(abs(open_last_section) > iqr(open_last_section) * 20,1);
open_end = length(aud_op) - open_end_frame + open_end;

length_open = open_end - open_start;
time_total = length_open/SF;
true_avi_SF = (avi_end - avi_start)/time_total;

SF_ratio = true_avi_SF/SF;
time_vector_avi = [1:length(aud_avi)]' - avi_start;
time_vector_avi = time_vector_avi + (floor(SF_ratio * open_start));







%% create output structures
Avisoft_times.avi_start = avi_start;
Avisoft_times.avi_end = avi_end + 1;
Avisoft_times.alignedFrames = time_vector_avi;
Avisoft_times.SF = true_avi_SF;
Avisoft_times.open_ephys_start = open_start;
Avisoft_times.open_ephys_end = open_end + 1;
Avisoft_times.sampling_ratio = SF_ratio;
end

