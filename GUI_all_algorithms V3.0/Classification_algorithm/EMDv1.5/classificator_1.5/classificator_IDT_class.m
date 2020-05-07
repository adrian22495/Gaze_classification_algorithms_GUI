% I-DT model classificator
classdef classificator_IDT_class <  eye_tracker_raw_data_reader_class & ...             % Reader from eye tracker data
                                    eye_records_class & ...                             % Basic class for placing eye tracker data
                                    eye_tracker_raw_data_converter_ETU_degree & ...     % Convertor between ETU and degress in data
                                    eye_tracker_raw_data_filter_class & ...             % Eye tracker data filtering by range of degrees
                                    classificator_merge_class & ...                     % Creates sequences of eye movements
                                    classificator_saccade_amplitude_filter_class & ...  % Filtered saccades based theire amplitude
                                    classificator_datafile_output_class & ...           % Output sequences to the files
                                    classificator_get_percentage_class & ...            % Calculate percentage of movements of every type
                                    classificator_enumerations_class & ...              % Basic enumerations definitions
                                    classificator_time_rate_class & ...                 % Time step and sample rate definitions
                                    handle
    % This is implementation of I-DT classification
    % For classification you have to do next:
    % 1. Setup delta_t_sec property calling set.delta_t_sec method
    % 2. Setup dispersion_duration_sec_threshold property calling
    % set.dispersion_duration_sec_threshold
    % 3. Setup dispersion_threshold property calling
    % set.dispersion_threshold
    % 4. Call classify method to make classification
    
    properties (Hidden)
% Publicly accessible properties
        dispersion_duration_sec_threshold;  % Length of window in sec.
        dispersion_threshold;               % Value of dispersion threshold of I-DT classification model
% Hidden properties for inner use
        window_length;                      % Length of sample window for I-DT model
    end

    methods
% Classification function
        function classify(obj)
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Begin data classification with I-DT classifier in :',datestr(now),'\n'));
            end
            obj.calculate_delta_t();
% Setup initial window parameters
            [IDT_window_start IDT_window_ended] = obj.get_window_parameters(1);
            while (IDT_window_start < length(obj.eye_records) )
                IDT_dispersion = obj.get_window_dispersion(IDT_window_start, IDT_window_ended);
                if( IDT_dispersion < obj.dispersion_threshold)
% According I-DT classification we should add records to the window until
% its dispersion became equal to threshold
                    while ( (IDT_dispersion < obj.dispersion_threshold) && (IDT_window_ended < length(obj.eye_records)) )
                        IDT_window_ended = IDT_window_ended + 1;
                        IDT_dispersion = obj.get_window_dispersion(IDT_window_start, IDT_window_ended);
                    end
% Now we should mark all records inside window as a fixations
                    obj.eye_records( IDT_window_start:IDT_window_ended,obj.MOV_TYPE ) = obj.FIXATION_TYPE;
% Clear the window
                    [IDT_window_start IDT_window_ended] = obj.get_window_parameters(IDT_window_ended+1);
                else
% According I-DT classification first record in window is a saccade. So we
% should mark this first as saccade
                    obj.eye_records( IDT_window_start,obj.MOV_TYPE ) = obj.SACCADE_TYPE;
% Due to our method of velocity calculations we should mark previous point
% as saccade.
                    if( IDT_window_start > 1)
                        obj.eye_records( IDT_window_start-1,obj.MOV_TYPE ) = obj.SACCADE_TYPE;
                    end
% Remove first point from current window                    
                    IDT_window_start = IDT_window_start + 1;
                end
            end
% Now mark all invalid data as noise
            obj.eye_records( (obj.eye_records(:,obj.VALIDITY) == obj.DATA_INVALID ),obj.MOV_TYPE ) = obj.NOISE_TYPE;
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Complete data classification with I-DT classifier in :',datestr(now),'\n'));
            end
        end

% Access interface to publicly accessible properties
        function set.dispersion_duration_sec_threshold(obj, value)
            obj.dispersion_duration_sec_threshold = value;
        end

        function set.dispersion_threshold(obj, value)
            obj.dispersion_threshold = value;
        end

% Access to class service routines
        function result = get.window_length(obj)
            result = ceil(obj.sample_rate * obj.dispersion_duration_sec_threshold);
        end

% Calculating begin and end of sample window for current position
        function [window_start window_ended] = get_window_parameters(obj, current_record)
            window_start = current_record;
            window_ended = window_start + obj.window_length;
            if( window_ended > length(obj.eye_records) )
                window_ended = length(obj.eye_records);
            end
        end

% Calculate dispersion of window with specified starting and ending
% position
        function result = get_window_dispersion(obj, window_start, window_end)
% Looking for maximal and minimal x and y coordinates of the points in the window
                IDT_max_x = max ( obj.eye_records( window_start:window_end, obj.X_COORD) );
                IDT_max_y = max ( obj.eye_records( window_start:window_end, obj.Y_COORD) );
                IDT_min_x = min ( obj.eye_records( window_start:window_end, obj.X_COORD) );
                IDT_min_y = min ( obj.eye_records( window_start:window_end, obj.Y_COORD) );
% Calculate the dispersion of the points in the window
                result = abs( IDT_max_x - IDT_min_x ) + abs( IDT_max_y - IDT_min_y );
        end
    end
    
end
