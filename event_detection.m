function [times] = event_detection(data,timeBefore)


%sigmaTimeSec = 0.08;
%sigmaSamples = round(sigmaTimeSec*30000);
%window = fspecial('gaussian',[1, sigmaSamples*6],sigmaSamples);
%data = conv(data,window,'same');
data = abs(data);

figure;
plot(data);
yline(iqr(data))
thr = input(['black line = ' num2str(iqr(data)),'; input threshold value: '])
close all
figure;
plot(data)
yline(thr);
thr = input('reenter threshold or correct')


above_thr = find(data > thr);
below_thr = find(data < thr)+1;
times = intersect(above_thr',below_thr');

timeBefore = timeBefore * 30000;
times((times - timeBefore < 1)) = [];
wind_mat = (ones(length(times),timeBefore) .* [-timeBefore:1:-1]) + times;
wind_dat_onset = abs(data(wind_mat));
[rows,~] = find(wind_dat_onset > (thr)); rows = unique(rows);
times(rows) = [];


to_plot = 0; % If you want to varify the detection on a plot make 1
switch to_plot
    case 0
    case 1
        figure; plot ((data),'LineWidth',0.5);
        hold on
        for i=1:10%length(whisk_onset)
            xline(times(i),'r');
        end
        g = ones(length(times(5:30)),60000);
        g = g.*(-29999:30000);
        gx = g + times(5:30);
        gd = data(gx);
        figure; hold on; plot(mean(gd))
        %x = 1;
end




end

