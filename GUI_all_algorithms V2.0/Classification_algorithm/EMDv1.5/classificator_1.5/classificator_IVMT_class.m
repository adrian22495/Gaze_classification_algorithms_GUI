% This code is based on Mahmoud Mechehoul's code
% Department of Computer Science, Texas State University
% mm2026@txstate.edu

classdef classificator_IVMT_class < eye_tracker_raw_data_reader_class & ...             % Reader from eye tracker data
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
    % This is implementation of I-VT classification
    % For classification you have to do next:
    % 1. Setup delta_t_sec property calling set.delta_t_sec method
    % 2. Setup saccade_threshold property calling set.saccade_threshold method
    % 3. Call classify method to make classification
    properties (Hidden)
% Publicly accessible properties
        saccade_detection_threshold;            % saccade detection threshold value - see I-VT model description for details
        tmp_window_length;                      % Length of window for San Agustin algorithm, samples
        san_agustin_threshold = 0.1;            % Threshold value for San Agustin Idea
        window_duration_sec_threshold = 0.05;   % Length of window length for San Agustin algorithm, seconds
    end    
    
    methods
% Classification functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
% ================== I-VT + San Agustin idea ===================== %
%                                                                  %
%==================================================================%
        function classify(obj)
            if( obj.debug_mode ~= 0), fprintf(strcat('Begin data classification with I-VMT classifier in :',datestr(now),'\n')); end;
% ===============================================
%    STAGE 1 - Saccade/Non saccade classification
% ===============================================
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

% Perform a first stage classification - mark all saccade samples
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
% Now we move towards the end os samples until we meet sample of different
% type or end of the samples
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
% Now we have to expand all survived saccades to the one point from left
            tmp_type = obj.eye_records(2:end,obj.MOV_TYPE);
            tmp_type(length(tmp_type)+1) = NaN;
            obj.eye_records( ((obj.eye_records(:,obj.VALIDITY) == obj.DATA_VALID) & (tmp_type(:) == obj.SACCADE_TYPE)),obj.MOV_TYPE) = obj.SACCADE_TYPE;
% =================================================
%                  END OF STAGE 1
% =================================================

% =================================================
%      STAGE 2 - San Agustin idea
% =================================================
% From this point we ready to perform second stage classification. For all
% samples marked as fixations we should check the distribution of angles
% inside of every temporal window.
%
% Initialize temporal window
            TMP_window_begin = 1;
% Until we reach the end of eye positions records
            while (TMP_window_begin < length(obj.eye_records) )
% Calculate the offset of temporal window
                TMP_window_end = min( TMP_window_begin + obj.tmp_window_length, length(obj.eye_records) );
% Clear angle array
                angle_array = [];
                cX=[];
                cY=[];
                distances = [];
% Perform evaluation of all angles
                for i=TMP_window_begin+1:TMP_window_end
% Estimate angle for line built for i-1:i points
                    if( obj.eye_records(i-1, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE && obj.eye_records(i, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE )
                        x1 = obj.eye_records(i-1, obj.X_COORD); y1 = obj.eye_records(i-1, obj.Y_COORD);
                        x2 = obj.eye_records(i, obj.X_COORD);   y2 = obj.eye_records(i, obj.Y_COORD);
                        dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
                        if( isempty( distances ) ), distances(1) = dist;
                        else distances(end+1)= dist;
                        end
                    end
                end
                
                for i=TMP_window_begin+1:TMP_window_end
% Estimate angle for line built for i-1:i points
                    if( obj.eye_records(i-1, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE && obj.eye_records(i, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE )
                        x1 = obj.eye_records(i-1, obj.X_COORD); y1 = obj.eye_records(i-1, obj.Y_COORD);
                        x2 = obj.eye_records(i, obj.X_COORD);   y2 = obj.eye_records(i, obj.Y_COORD);
                        dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
                        %if( dist >= mean(distances) )
                            if( x1==x2 ), current_angle = pi/2; if( y1>y2), current_angle=3*pi/2; end
                            else current_angle = atan( (y2-y1)/(x2-x1) ); if( x1>x2 ), current_angle=current_angle+pi; end
                            end
                            if(current_angle>2*pi ), current_angle=current_angle-2*pi;end;
                            if( isempty( angle_array ) ), angle_array(1) = current_angle;
                            else angle_array(end+1)=current_angle;
                            end
% Change angles into cartesian circle R=1
                            ccx = sin(current_angle);
                            ccy = cos(current_angle);
                            if( isempty( cX ) ), cX(1) = ccx; else cX(end+1) = ccx; end;
                            if( isempty( cY ) ), cY(1) = ccy; else cY(end+1) = ccy; end;
                        %end
                    end
                end


% Remove all repititive values from array of angles
                %angle_array = angle_array + pi/2;
                %uniq_angles = unique(angle_array);
                %for i=2:length(angle_array), changes_array(i) = angle_array(i) - angle_array(i-1); end
% Move to next temporal window
                if( sqrt( mean(cX)^2+mean(cY)^2) > obj.san_agustin_threshold )
                    for i=TMP_window_begin+1:TMP_window_end
                        if (obj.eye_records(i, obj.MOV_TYPE ) ~= obj.SACCADE_TYPE ), obj.eye_records(i, obj.MOV_TYPE ) = obj.PURSUIT_TYPE; end;
                    end
                end
                
                TMP_window_begin = TMP_window_end + 1;
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
% find its end and check if the next sample after end if fixation or
% saccade
                if( obj.eye_records( current_position, obj.MOV_TYPE ) == obj.PURSUIT_TYPE )
% First search for end of pursuit
                    while( current_position < length( obj.eye_records) && obj.eye_records( current_position, obj.MOV_TYPE ) == obj.PURSUIT_TYPE )
                        current_position = current_position + 1;
                    end
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

            if( obj.debug_mode ~= 0), fprintf(strcat('Complete data classification with I-VMT classifier in :',datestr(now),'\n')); end
        end
        
% Access interface to publicly accessible properties
        function set.saccade_detection_threshold(obj,value), obj.saccade_detection_threshold = value; end
        function set.window_duration_sec_threshold(obj, value), obj.window_duration_sec_threshold = value; end
        function result = get.tmp_window_length(obj), result = ceil(obj.sample_rate * obj.window_duration_sec_threshold); end
        function set.san_agustin_threshold(obj,value), obj.san_agustin_threshold = value; end
       
    end

end
