classdef auditory_ephys_exp
    % Processed data from ephys experiment for auditory-tectile project
    % starting Dec 2020
    %   This class contains the unit data together with the relevent analog
    %   channels and timestamps of digital events. 
    %   The class contains the following info:
    %       Units firing times in frames
    %       Units identity - good,mua
    %       Start and end time for experiment (From avi synchronization)
    %       Avisoft timestamps
    %       Cameras timestamps
    %       Experiment conditions times
    %       Encoder data
    %       Micropone channel (either from avi or from openephys)
    %       Speaker activity times
    %       Whisking signal (on/off)
    %       DLC whisking signal
    
    properties
        experiment_ID
        Units   % Units times and identity !!! align with avi if needed
        Mic     % Microphone signal from avisoft (if exist)
        Whisking_simple     % Whisking on/off signal 
        Whisking_DLC    % Complex whisking signal etected with DLC
        Encoder     % Running data of mouse
        Speaker     % Times of evoked speaker signal
        Conditions  % Times of different experiment conditions
        Avisoft_times   % Start, end, SF, and timestamps of mic signal from Avisoft
        Cams    % Timestamps of camrea frames !!!align with avi if needed
        SF % experiment sampling rate
        firing_per_condition
        alum_prediction
       
    end
    
    properties (Hidden = true)
        parameters_file % the full path for the parameters file
        color % The RGB colors for the figures
        experiment_path % The full path to the experiment folder
    end
    
    methods 
        function obj = auditory_ephys_exp(parameters_file) 
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            run(parameters_file);
            obj.experiment_ID = expName;
            obj.experiment_path = exp_path;
            obj.parameters_file = parameters_file;
            obj.SF = SF;
            cd(exp_path);
            try
                whisk_filename = dir('*whisking.csv');
                obj.Whisking_simple.Raw = csvread([whisk_filename.folder,'/',whisk_filename.name]);
            catch
            end
            obj.color.Aluminum = [0.0039 * 210 , 0.0039 * 30 , 0.0039 * 75]; % D21E4B adobe color wheel
            obj.color.noObject = [0.0039 * 232 , 0.0039 * 134 , 0.0039 * 65]; % 284B59 adobe color wheel
            obj.color.Muted = [0.0039 * 40, 0.0039 * 75, 0.0039 * 89]; % E88641 adobe color wheel
            obj.experiment_path = exp_path;
            obj.Cams.Whisking = syncCamera(obj.parameters_file,1);
            obj.Cams.Face = syncCamera(obj.parameters_file,2);
            obj.Units = sorted_from_phy(obj.parameters_file);
            %% get the data from openephys
            data_temp = load_open_ephys_binary([exp_path,'/experiment1/recording1/structure.oebin'],'continuous',1,'mmap');
            obj.Encoder = double(data_temp.Data.Data.mapped(32 + encoder,:));
            try
                [obj.Avisoft_times] = avisoftSync(obj.parameters_file);
                cd(exp_path);
                files = dir('*.wav');
                obj.Mic = audioread([files.folder,'/',files.name]);
            catch
                warning("Can't run avisoft sync, assigning MIC signal from open ephys");
                obj.Mic = double(data_temp.Data.Data.mapped(32 + mic,:));
                obj.Avisoft_times.open_ephys_end = length(obj.Mic);
                
            end
            obj.Conditions = obj.create_condition_times;
            %obj.firing_per_condition = obj.firing_rates_basic;
            
        end
        
        function firing_rates_mean = firing_rates_basic(obj)
            %Creates a matrix with all units firing rate average for each
            %condition normalized by the amount of whisking
            
            % reshape unit matrices into a 3d array 
            good_units = obj.Units.good.times;
            good_units = reshape(good_units,[1,size(good_units,2),size(good_units,1)]);
            mua_units = obj.Units.mua.times;
            mua_units = reshape(mua_units,[1,size(mua_units,2),size(mua_units,1)]);
            
            
            % find all firing events for each condition and sum it
            good_units_conditions = good_units > obj.Conditions.start_end_times(:,1) & good_units < obj.Conditions.start_end_times(:,2);
            good_units_conditions = sum(good_units_conditions,2);
            good_units_conditions = reshape(good_units_conditions,[size(good_units_conditions,1),size(good_units_conditions,3)]);
            mua_units_conditions = mua_units > obj.Conditions.start_end_times(:,1) & mua_units < obj.Conditions.start_end_times(:,2);
            mua_units_conditions = sum(mua_units_conditions,2);
            mua_units_conditions = reshape(mua_units_conditions,[size(mua_units_conditions,1),size(mua_units_conditions,3)]);
            
            % add firing events for instances of the condition and
            % normalize by the total time of each condition
            total_times_conditions = obj.Conditions.start_end_times(:,2) - obj.Conditions.start_end_times(:,1);
            times_conditions = accumarray(obj.Conditions.condition_classification,total_times_conditions)/obj.SF;
            mua = zeros(size(mua_units_conditions,2),length(times_conditions))';
            good = zeros(size(good_units_conditions,2),length(times_conditions))';
            for unit = 1:size(mua,2)
                mua(:,unit) = accumarray(obj.Conditions.condition_classification,mua_units_conditions(:,unit));
            end
            mua = mua./times_conditions;
            for unit = 1:size(good,2)
                good(:,unit) = accumarray(obj.Conditions.condition_classification,good_units_conditions(:,unit));
            end
            good = good./times_conditions;
            
            %% extract total whisking time in time windows 
            % find the equivelent index for the frames of each condition
            [d(:,1),ix(:,1)] = min(abs(obj.Cams.Whisking.csv_aligned_frames - obj.Conditions.start_end_times(:,1)'));
            [d(:,2),ix(:,2)] = min(abs(obj.Cams.Whisking.csv_aligned_frames - obj.Conditions.start_end_times(:,2)'));
            logical_whisking = logical(obj.Whisking_simple.Raw);
            for condition=1:length(ix)
                total_whisking(1,condition) = sum(logical_whisking(ix(condition,1):ix(condition,2)));
            end
            total_whisking = total_whisking';
            total_time_whisking = accumarray(obj.Conditions.condition_classification,total_whisking)/obj.Cams.Whisking.fps;
            mua_normalized = mua./total_time_whisking;
            good_normalized = good./total_time_whisking;
            
            firing_rates_mean.mua = mua_normalized;
            firing_rates_mean.good = good_normalized;
            firing_rates_mean.Whisking_per_condition = total_time_whisking;
            firing_rates_mean.Time_in_condition = times_conditions;
        end
        
        
        function clean_condition_start_end = create_condition_times(obj)
            %Extracts all the spikes in each condition
            % Run loop for the three conditions and find all spikes between
            % times
            load([obj.experiment_path,'/condition_extract.mat']);
            clean_condition_start_end.condition_classification = condition_extract.condition_classification;
            clean_condition_start_end.start_end_times(:,1) = double(condition_extract.all_changes);
            clean_condition_start_end.start_end_times(1:end-1,2) = double(condition_extract.all_changes(2:end));
            clean_condition_start_end.start_end_times(end,2) = obj.Avisoft_times.open_ephys_end;
            clean_condition_start_end.start_end_times(:,1) = clean_condition_start_end.start_end_times(:,1) + (condition_extract.artifact_times(2)/1000 * obj.SF);
            clean_condition_start_end.start_end_times(:,2) = clean_condition_start_end.start_end_times(:,2) - (condition_extract.artifact_times(1)/1000 * obj.SF);
            clean_condition_start_end.artifact_times = condition_extract.artifact_times;
            if clean_condition_start_end.start_end_times(end,1) > clean_condition_start_end.start_end_times(end,2)
               clean_condition_start_end.start_end_times(end,:) = [];
               clean_condition_start_end.condition_classification(end) = [];
            end
        end
        
        function prediction = createTableForPrediction(obj,model)
            %% load all relevant information
            cd('/home/ben/Z')
            load(['All_whiskers_videos/Analyzed/',obj.experiment_ID,'.mat'],'motion_1','motion_2','motSVD_1','motSVD_2');

            %% 
            timeAfter = obj.Cams.Whisking.fps * 0.05-1;
            totalFrameWindow = 1:obj.Cams.Whisking.fps * 0.06;
            timeMatrix = (obj.Cams.Whisking.csv_start_frame:length(obj.Cams.Whisking.csv_aligned_frames) - timeAfter)';
            timeMatrix = totalFrameWindow + timeMatrix-1;
        %% create the features for model
            predictionTable = table;
            motion_1_labeled = motion_1(timeMatrix);
            motion_2_labeled = motion_2(timeMatrix);
            predictionTable.motion_1_mean = mean(motion_1_labeled,2);
            predictionTable.motion_2_mean = mean(motion_2_labeled,2);
            predictionTable.motion_1_STD = std(motion_1_labeled,0,2);
            predictionTable.motion_2_STD = std(motion_2_labeled,0,2);
            if obj.Cams.Whisking.fps == 300
                motion_1_labeled = single(resample(double(motion_1_labeled),400,300,'Dimension',2));
                motion_2_labeled = single(resample(double(motion_2_labeled),400,300,'Dimension',2));
            end
            predictionTable.motion_1 = motion_1_labeled;
            predictionTable.motion_2 = motion_2_labeled;

            for mot = 1:10
                tempMot_1 = motSVD_1(:,mot);
                tempMot_2 = motSVD_2(:,mot);
                tempMot_1_resample = tempMot_1(timeMatrix);
                tempMot_2_resample = tempMot_2(timeMatrix);
                if obj.Cams.Whisking.fps == 300
                    tempMot_1_resample = single(resample(double(tempMot_1_resample),400,300,'Dimension',2));
                    tempMot_2_resample = single(resample(double(tempMot_2_resample),400,300,'Dimension',2));     
                end
    
                predictionTable.(['MotSVD_1_',num2str(mot)]) = tempMot_1_resample;
                predictionTable.(['MotSVD_2_',num2str(mot)]) = tempMot_2_resample;
            end
            zeroTime = timeMatrix(:,5);
            zeroTime(:,2) = obj.Cams.Whisking.csv_aligned_frames(zeroTime);
            yfit = model.predictFcn(predictionTable);
            alumFit = zeroTime(yfit == 1,:);
            tol = 3000/max(abs(alumFit(:,2)));
            prediction.all = uniquetol(alumFit(:,2),tol);
            
            %% get predictions for each condition
            alumTimes = obj.Conditions.start_end_times(obj.Conditions.condition_classification == 1,:);
            mutedTimes = obj.Conditions.start_end_times(obj.Conditions.condition_classification == 2,:);
            nanTimes = obj.Conditions.start_end_times(obj.Conditions.condition_classification == 3,:);
            
            prediction.Aluminum = [];
            prediction.Muted = [];
            prediction.noObject = [];
            
            
            for i = 1:length(alumTimes)
                tempPred = prediction.all(prediction.all > alumTimes(i,1) & prediction.all < alumTimes(i,2));
                prediction.Aluminum = [prediction.Aluminum;tempPred];
            end
            
            for i = 1:length(mutedTimes)
                tempPred = prediction.all(prediction.all > mutedTimes(i,1) & prediction.all < mutedTimes(i,2));
                prediction.Muted = [prediction.Muted;tempPred];
            end
            
            for i = 1:length(nanTimes)
                tempPred = prediction.all(prediction.all > nanTimes(i,1) & prediction.all < nanTimes(i,2));
                prediction.noObject = [prediction.noObject;tempPred];
            end
            
        end


        
        
        
    end
end

