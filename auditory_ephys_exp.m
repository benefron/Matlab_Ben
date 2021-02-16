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
        OpenEpyhs_times     % Start, end and timestamps of all signals from openephys
        Cams    % Timestamps of camrea frames !!!align with avi if needed
       
    end
    
    properties (Hidden = true)
        parameters_file % the full path for the parameters file
        color % The RGB colors for the figures
        experiment_path % The full path to the experiment folder
    end
    
    methods 
        function obj = create_experiment(parameters_file) 
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            run(parameters_file);
            obj.experiment_ID = expName;
            obj.parameters_file = parameters_file;
            obj.color.Aluminum = [0.0039 * 210 , 0.0039 * 30 , 0.0039 * 75]; % D21E4B adobe color wheel
            obj.color.noObject = [0.0039 * 232 , 0.0039 * 134 , 0.0039 * 65]; % 284B59 adobe color wheel
            obj.color.Muted = [0.0039 * 40, 0.0039 * 75, 0.0039 * 89]; % E88641 adobe color wheel
            obj.experiment_path = exp_path;
            %% get the synced times and vectors from avi and camera
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

