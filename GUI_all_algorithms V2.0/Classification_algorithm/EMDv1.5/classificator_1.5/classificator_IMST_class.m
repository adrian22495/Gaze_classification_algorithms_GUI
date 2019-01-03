classdef classificator_IMST_class < eye_tracker_raw_data_reader_class & ...             % Reader from eye tracker data
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
    % This is implementation of I-MST classification
    % For classification you have to do next:
    % 1. Setup delta_t_sec property calling set.delta_t_sec method
    % 2. Setup saccade_detection_threshold property calling set.saccade_detection_threshold method
    % 3. Setup window_size property calling set.window_size method
    % 4. Call classify method to make classification
    
    properties (Hidden)
% Publicly accessible properties
        saccade_detection_threshold;    % Saccade detection threshold value - see I-VT model description for details
        window_size;                    % Size of one batch in elements;
% Hidden properties for inner use
        window_count;                   % Count of batches in current data set
    end
    
    methods
% Classification function
        function classify(obj)
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Begin data classification with I-MST classifier in :',datestr(now),'\n'));
            end
            obj.calculate_delta_t();
% Begin processing data records by moving windows
            for k=1:obj.window_count
% Get begin, end, and length of the current window
                [start_record ended_record window_length] = obj.get_window_parameters(k);
% Computing Euclidean distance between each eye position inside the window
                X = zeros(window_length,2);
                X(:,1) = obj.eye_records(start_record:ended_record, obj.X_COORD);
                X(:,2) = obj.eye_records(start_record:ended_record, obj.Y_COORD);
                distance_matrix = squareform(pdist(X,'euclidean'));
                distance_matrix(logical(eye(size(distance_matrix)))) = inf;

% Building minimum spanning tree using Prim's algorithm. I'm afraid that
% this implementation is not acceptable for large data arrays. It's better
% to use something to improve this.
                vertices = zeros( window_length,1 );    % List of visited vertices
                vertices(1) = 1;                        % Mark first vertex as visited
                count_v = 1;                            % Total amount of visited vertices
                while (count_v<window_length)           % While some of the vertices is unvisited
                    min_distance = inf;                 % Searching for minimal distance betwenn visited and unvisited vertices
                    for i=1:window_length
                        if (vertices(i) == 1)
                            for j=1:window_length
                                if (vertices(j) == 0)
                                    if( min_distance > distance_matrix( i,j ) )
                                        min_distance = distance_matrix( i,j );
                                        v1 = i;
                                        v2 = j;
                                    end
                                end
                            end
                        end
                    end
                    vertices(v2) = 1;                   % Mark selected unvisited vertex as visited
% Comparing new edge of the MST with threshold value
                    if( min_distance < obj.saccade_detection_threshold )
% If length of this edge is less than threshold value then fixation
% detected but only for valid data
                        if( obj.eye_records( v1 + start_record -1,obj.VALIDITY) == obj.DATA_VALID)
                            obj.eye_records( v1 + start_record - 1 ,obj.MOV_TYPE ) = obj.FIXATION_TYPE;
                        else
                            obj.eye_records( v1 + start_record - 1 ,obj.MOV_TYPE ) = obj.NOISE_TYPE;
                        end
                        if( obj.eye_records( v2 + start_record -1,obj.VALIDITY) == obj.DATA_VALID)
                            obj.eye_records( v2 + start_record - 1 ,obj.MOV_TYPE ) = obj.FIXATION_TYPE;
                        else
                            obj.eye_records( v2 + start_record - 1 ,obj.MOV_TYPE ) = obj.NOISE_TYPE;
                        end
                    else
% saccade detected but only for valid data
                        if( obj.eye_records( v1 + start_record -1,obj.VALIDITY) == obj.DATA_VALID)
                            obj.eye_records( v1 + start_record - 1 ,obj.MOV_TYPE ) = obj.SACCADE_TYPE;
                        else
                            obj.eye_records( v1 + start_record - 1 ,obj.MOV_TYPE ) = obj.NOISE_TYPE;
                        end
                        if( obj.eye_records( v2 + start_record -1,obj.VALIDITY) == obj.DATA_VALID)
                            obj.eye_records( v2 + start_record - 1 ,obj.MOV_TYPE ) = obj.SACCADE_TYPE;
                        else
                            obj.eye_records( v2 + start_record - 1 ,obj.MOV_TYPE ) = obj.NOISE_TYPE;
                        end
                    end
                    count_v = count_v + 1;
                end
            end
% Due to our method of velocity computations (left difference) we have to
% add left edge to our saccades
            tmp_type = obj.eye_records(2:end,obj.MOV_TYPE);
            tmp_type(length(tmp_type)+1) = NaN;
            obj.eye_records( ((obj.eye_records(:,obj.VALIDITY) == obj.DATA_VALID) & (tmp_type(:) == obj.SACCADE_TYPE)),obj.MOV_TYPE) = obj.SACCADE_TYPE;

            if( obj.debug_mode ~= 0)
                fprintf(strcat('Complete data classification with I-MST classifier in :',datestr(now),'\n'));
            end
        end

% Public access interface for class properties
        function set.saccade_detection_threshold(obj,value)
            obj.saccade_detection_threshold = value;
        end

        function set.window_size(obj, value)
            obj.window_size = value;
        end

        function result = get.window_count(obj)
            result = ceil(length (obj.eye_records) / obj.window_size );
        end

        function [start_window end_window window_length] = get_window_parameters(obj, count)
            start_window = obj.window_size * (count-1) + 1;
            end_window = start_window + obj.window_size - 1;
            if (end_window > length(obj.eye_records) )
                end_window = length( obj.eye_records);
            end
            window_length = end_window - start_window + 1;
        end
    end

end
