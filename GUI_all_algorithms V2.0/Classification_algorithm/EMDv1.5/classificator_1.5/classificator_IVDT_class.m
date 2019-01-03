% This code is based on Mahmoud Mechehoul's code
% Department of Computer Science, Texas State University
% mm2026@txstate.edu

classdef classificator_IVDT_class < eye_tracker_raw_data_reader_class & ...             % Reader from eye tracker data
                                    eye_records_class & ...                             % Basic class for placing eye tracker data
                                    eye_tracker_raw_data_converter_ETU_degree & ...     % Convertor between ETU and degress in data
                                    eye_tracker_raw_data_filter_class & ...             % Eye tracker data filtering by range of degrees
                                    classificator_merge_class & ...                     % Creates sequences of eye movements
                                    classificator_saccade_amplitude_filter_class & ...  % Filtered saccades based theire amplitude
                                    classificator_datafile_output_class & ...           % Output sequences to the files
                                    classificator_get_percentage_class & ...            % Calculate percentage of movements of every type
                                    classificator_enumerations_class & ...              % Basic enumerations definitions
                                    classificator_time_rate_class & ...                 % Time step and sample rate definition
                                    handle
    % This is implementation of I-VDT classification
    % For classification you have to do next:
    % 1. Setup delta_t_sec property calling set.delta_t_sec method
    % 2. Setup dispersion_duration_sec_threshold property calling set.idt_dispersion_threshold
    % 3. Setup dispersion_threshold property calling set.idt_window_length
    % 4. Call classify method to make classification
    properties (Hidden)
% Publicly accessible properties
        saccade_detection_threshold;            % saccade detection threshold value, deg/sec - see I-VT model description for details
        dispersion_duration_sec_threshold = 0.1;% Length of window in sec.
        idt_dispersion_threshold = 1.0;         % Threshold value for maximal allowed dispersion of fixation, degrees
        idt_window_length;                      % Length of sample window for I-DT model
    end    
    
    methods
% Classification functions
        function classify(obj)
            if( obj.debug_mode ~= 0), fprintf(strcat('Begin data classification with I-VDT classifier in :',datestr(now),'\n')); end
% =================================================================
%    STAGE 1 - Saccade/non saccade classification, I-VT algorithm
% =================================================================
            obj.calculate_delta_t();
            x_velocity_degree = zeros( length(obj.eye_records),1 );
            y_velocity_degree = zeros( length(obj.eye_records),1 );
% Calculate absolute degree velocity of our records
            x_velocity_degree( 2:end ) =(   obj.eye_records( 2:end,obj.X_COORD ) - ...
                                            obj.eye_records( 1:end-1,obj.X_COORD ) ) / obj.delta_t_sec;
            y_velocity_degree( 2:end ) =(   obj.eye_records( 2:end,obj.Y_COORD ) - ...
                                            obj.eye_records( 1:end-1,obj.Y_COORD ) ) / obj.delta_t_sec;
% First point is a special case
            x_velocity_degree(1) = 0;
            y_velocity_degree(1) = 0;
            obj.eye_records(:,obj.VELOCITY) = sqrt( x_velocity_degree.^2 + y_velocity_degree.^2 );
% First point is a special case
            obj.eye_records(1,obj.MOV_TYPE ) = obj.NOISE_TYPE;
            obj.eye_records(1,obj.VELOCITY) = 0;

% Perform a first stage classification - mark all saccade's samples
            obj.eye_records( ((abs(obj.eye_records(:,obj.VELOCITY)) >= obj.saccade_detection_threshold) & (obj.eye_records(:,obj.MOV_TYPE) ~= obj.NOISE_TYPE)),obj.MOV_TYPE ) = obj.SACCADE_TYPE;
% And all other samples as fixations
            obj.eye_records( ((abs(obj.eye_records(:,obj.VELOCITY)) < obj.saccade_detection_threshold) & (obj.eye_records(:,obj.MOV_TYPE) ~= obj.NOISE_TYPE)),obj.MOV_TYPE ) = obj.FIXATION_TYPE;
% =================================================
%   STAGE 1 - Prelimenary filtration of saccades
% =================================================
% Perform a preliminary filtration - search all saccades, check their
% amplitude and length and mark them back as nonclassified if failed.
            i=1;
            while( i < length(obj.eye_records) )
% Get the first position of our classified eye movement
                onset_position = i;
% Now we move towards the end of samples sequence until we meet sample of
% different type or end of the samples
                while ( i<length(obj.eye_records) && obj.eye_records( i , obj.MOV_TYPE) == obj.eye_records( onset_position , obj.MOV_TYPE )), i=i+1; end
% Now we found our offset position
                offset_position = i;
% Now we should exclude saccades with amplitudes less then 0.5 degree
% And saccades that outside of allowed range, if necessary
                if ( obj.eye_records( onset_position, obj.MOV_TYPE ) == obj.SACCADE_TYPE )
                    saccade_amplitude_x =   max( obj.eye_records( onset_position:offset_position,obj.X_COORD ) ) - ...
                                            min( obj.eye_records( onset_position:offset_position,obj.X_COORD ) );
                    saccade_amplitude_y =   max( obj.eye_records( onset_position:offset_position,obj.Y_COORD ) ) - ...
                                            min( obj.eye_records( onset_position:offset_position,obj.Y_COORD ) );                                        
                    saccade_amplitude = sqrt ( saccade_amplitude_x^2 + saccade_amplitude_y^2 );
                    saccade_length = (offset_position-onset_position+1) ;

                    if( saccade_amplitude < obj.minimal_saccade_amplitude || ...
                        saccade_length < obj.minimal_saccade_length )
                        obj.eye_records( onset_position:offset_position, obj.MOV_TYPE ) = obj.FIXATION_TYPE;
                    end
                end
            end
% Now we have to expand all non-filtered saccades to the one point from left
            tmp_type = obj.eye_records(2:end,obj.MOV_TYPE);
            tmp_type(length(tmp_type)+1) = NaN;
            obj.eye_records( ((obj.eye_records(:,obj.VALIDITY) == obj.DATA_VALID) & (tmp_type(:) == obj.SACCADE_TYPE)),obj.MOV_TYPE) = obj.SACCADE_TYPE;
% =================================================
%                  END OF STAGE 1
% =================================================

% =================================================
%      STAGE 2 - Fixation/Pursuit separation
%                I-DT algorithm
% =================================================
% Now we ready to perform second stage classification. All
% samples that were marked as fixation should be checked in order to determine
% if this is really fixation or smooth pursuit. For this purpose we are using
% I-DT classifier. It's an I-DT algorithm with one feature - we should left
% intact samples that was classified before us.
%
% Initialize I-DT window
            IDT_window_begin = 1;
            IDT_window_end = min( IDT_window_begin + obj.idt_window_length, length(obj.eye_records) );
% Until we reach the end of array
            while (IDT_window_begin < length(obj.eye_records) )
% Calculate dispersion  for current window
                IDT_dispersion = obj.get_idt_dispersion( IDT_window_begin, IDT_window_end);
% Check if dispersion of current window is lesser than threshold value
                if( IDT_dispersion <  obj.idt_dispersion_threshold )
% If true we begin increase the window by adding samples to it until it's
% dispersion becames equal to threshold value
                    while ( (IDT_dispersion < obj.idt_dispersion_threshold ) && (IDT_window_end < length(obj.eye_records)) )
                        IDT_window_end = IDT_window_end + 1;
                        IDT_dispersion = obj.get_idt_dispersion(IDT_window_begin, IDT_window_end);
                    end
% Now we successfully discover our fixation and should mark all samples as
% fixation
                    for i = IDT_window_begin : IDT_window_end
                        if( obj.eye_records( i, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE ), obj.eye_records( i, obj.MOV_TYPE ) = obj.FIXATION_TYPE; end
                    end
% And clean the window
                    IDT_window_begin = IDT_window_end+1;
                    IDT_window_end = min( IDT_window_begin + obj.idt_window_length, length(obj.eye_records) );
                else
% Otherwise (in case when dispersion of our window is larger than dispersion
% threshold value)
                    if( obj.eye_records( IDT_window_begin, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE ), obj.eye_records( IDT_window_begin, obj.MOV_TYPE ) = obj.PURSUIT_TYPE; end
% Mark previous point as pursuit in case if it wasn't saccade
                    if( IDT_window_begin > 1), if( obj.eye_records( IDT_window_begin-1, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE ), obj.eye_records( IDT_window_begin-1, obj.MOV_TYPE ) = obj.PURSUIT_TYPE; end; end
% And finally remove first sample from the window
                    IDT_window_begin = IDT_window_begin+1;
                end
            end
% =================================================
%                  END OF STAGE 2
% =================================================

% =================================================
%             STAGE 3 - Merge of pursuit
% =================================================
            current_position = 1;
            while( current_position < length(obj.eye_records) )
% Check if current sample is pursuit. If it really is pursuit we should
% find its end and check if the next sample after end is fixation or
% saccade
                if( obj.eye_records( current_position, obj.MOV_TYPE ) == obj.PURSUIT_TYPE )
% First search for end of pursuit
                    while( current_position < length( obj.eye_records) && obj.eye_records( current_position, obj.MOV_TYPE ) == obj.PURSUIT_TYPE ), current_position = current_position + 1; end
% Now we have to check what type is the next sample range
                    begin_position = current_position;
                    while(  current_position < length( obj.eye_records) && ...
                            obj.eye_records( current_position, obj.MOV_TYPE ) ~= obj.PURSUIT_TYPE &&...
                            obj.eye_records( current_position, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE )
                        current_position = current_position + 1;
                    end
                    if( current_position <= length( obj.eye_records) && ...
                        obj.eye_records( current_position, obj.MOV_TYPE ) == obj.PURSUIT_TYPE && ...
                        (current_position - begin_position) * obj.delta_t_sec * 1000 < 100 )
                        obj.eye_records( begin_position:current_position, obj.MOV_TYPE ) = obj.PURSUIT_TYPE;
                    end
                end
                current_position = current_position + 1;
            end
% =================================================
%                  END OF STAGE 3
% =================================================
% Now we have to expand all classified pursuits to the one point from left
            tmp_type = obj.eye_records(2:end,obj.MOV_TYPE);
            tmp_type(length(tmp_type)+1) = NaN;
            obj.eye_records( ((obj.eye_records(:,obj.VALIDITY) == obj.DATA_VALID) & (tmp_type(:) == obj.PURSUIT_TYPE)),obj.MOV_TYPE) = obj.PURSUIT_TYPE;

% Now we mark every invalid point as noise
            obj.eye_records( (obj.eye_records(:,obj.VALIDITY) == obj.DATA_INVALID),obj.MOV_TYPE ) = obj.NOISE_TYPE;
% And we marl every noise points as invalid
            obj.eye_records( (obj.eye_records(:,obj.MOV_TYPE) == obj.NOISE_TYPE),obj.VALIDITY ) = obj.DATA_INVALID;

% This is special case
            obj.eye_records( 1 , obj.MOV_TYPE ) = obj.eye_records( 2 , obj.MOV_TYPE );

            if( obj.debug_mode ~= 0), fprintf(strcat('Complete data classification with I-VDT classifier in :',datestr(now),'\n')); end
        end

% Access interface to publicly accessible properties
        function set.saccade_detection_threshold(obj,value), obj.saccade_detection_threshold = value; end
        function result = get.idt_window_length(obj), result = ceil(obj.sample_rate * obj.dispersion_duration_sec_threshold); end
        function set.dispersion_duration_sec_threshold(obj,value), obj.dispersion_duration_sec_threshold = value; end
        function set.idt_dispersion_threshold(obj,value), obj.idt_dispersion_threshold = value; end

% Get the dispersion for the given window
        function result = get_idt_dispersion( obj, onset_position, offset_position )
% Looking for maximal and minimal x and y coordinates of the points in the window
            IDT_max_x = max ( obj.eye_records( onset_position:offset_position, obj.X_COORD) );
            IDT_max_y = max ( obj.eye_records( onset_position:offset_position, obj.Y_COORD) );
            IDT_min_x = min ( obj.eye_records( onset_position:offset_position, obj.X_COORD) );
            IDT_min_y = min ( obj.eye_records( onset_position:offset_position, obj.Y_COORD) );
% Calculate the dispersion of the points in the window
            result = abs( IDT_max_x - IDT_min_x ) + abs( IDT_max_y - IDT_min_y );
        end

    end

end
