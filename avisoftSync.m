
function [Avisoft_times,aud_op] = avisoftSync(params_file)
%% load the files and create the variables
run(params_file);
data = load_open_ephys_binary([exp_path,'/experiment1/recording1/structure.oebin'],'continuous',1,'mmap');
channel_number = data.Header.num_channels;
aud_op = data.Data.Data.mapped(channel_number-(8-mic),:);
cd(exp_path);
files = dir('*.wav');
[aud_avi] = audioread([files.folder,'/',files.name]);

%% Find start and end point in both variables
aud_op = normalize(double(aud_op),"range",[-1,1]);

avi_thr = iqr(aud_avi) * 100;
op_thr =  iqr(aud_op) * 100;

avi_cross = find(aud_avi > avi_thr);
op_cross = find(aud_op > op_thr);



length_open = op_cross(end) - op_cross(1);
time_total = length_open/SF;
true_avi_SF = (avi_cross(end) - avi_cross(1))/time_total;

SF_ratio = true_avi_SF/SF;
time_vector_avi = [1:length(aud_avi)]' - avi_cross(1);
time_vector_avi = time_vector_avi + (floor(SF_ratio * op_cross(1)));







%% create output structures
Avisoft_times.avi_start = avi_cross(1);
Avisoft_times.avi_end = avi_cross(end) + 1;
Avisoft_times.alignedFrames = time_vector_avi;
Avisoft_times.SF = true_avi_SF;
Avisoft_times.open_ephys_start = op_cross(1);
Avisoft_times.open_ephys_end = op_cross(end) + 1;
Avisoft_times.sampling_ratio = SF_ratio;
end

