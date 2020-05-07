% This code is based on example from wikipedia
% http://en.wikipedia.org/wiki/Viterbi_algorithm and
% using part of David Conger's code
%
% Baum-Welch algorithm is based on S. Jayarathna, D.H Koh, and S.M. Gowda
% code.
% Department of Computer Science, Texas State University - San Marcos
% San Marcos
% sampath@txstate.edu, dk1132@txstate.edu, sm1449@txstate.edu

classdef classificator_IHMM_class < eye_tracker_raw_data_reader_class & ...             % Reader from eye tracker data
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
% This is implementation of I-HMM classification
    % For classification you have to do next:
    % 1. Setup delta_t_sec property calling set.delta_t_sec method
    % 2. Setup saccade_threshold property calling set.saccade_threshold method
    % 3. Setup viterbi_sample_size property calling set.viterbi_sample_size method
    % 4. Setup baum_welch_reiterations_count property calling set.baum_welch_reiterations_count method
    % 5. Call classify method to make classification
    
    properties (Hidden)
% Publicly accessible properties
        saccade_detection_threshold;    % Treshold value for initial classification
        viterbi_sample_size;            % Size of the window for Viterbi sampler
        baum_welch_reiterations_count;  % Count of iterations for Baum-Welch algorithm
% Hidden propertiesfor inner use
        window_count;                   % Count of batches in current data set
    end
    
    methods
% Classification function
        function classify(obj)
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Begin data classification with I-HMM classifier in :',datestr(now),'\n'));
            end
            obj.calculate_delta_t();

            x_velocity_degree = zeros( length(obj.eye_records),1 );
            y_velocity_degree = zeros( length(obj.eye_records),1 );
% Calculate absolute degree velocity of our records
            x_velocity_degree( 2:end ) =(   obj.eye_records( 2:end,obj.X_COORD ) - ...
                                            obj.eye_records( 1:end-1,obj.X_COORD) ) / obj.delta_t_sec;
            y_velocity_degree( 2:end ) =(   obj.eye_records( 2:end,obj.Y_COORD ) - ...
                                            obj.eye_records( 1:end-1,obj.Y_COORD) ) / obj.delta_t_sec;
% First point is a special case
            x_velocity_degree(1) = 0;
            y_velocity_degree(1) = 0;
            obj.eye_records(:,obj.VELOCITY) = sqrt( x_velocity_degree.^2 + y_velocity_degree.^2 );            
% First point is a special case
            obj.eye_records(1,obj.VELOCITY) = 0;
            
% Now we mark fixations
            obj.eye_records( (abs(obj.eye_records(:,obj.VELOCITY)) < obj.saccade_detection_threshold),obj.MOV_TYPE ) = obj.FIXATION_TYPE;
% Now we mark saccades
            obj.eye_records( (abs(obj.eye_records(:,obj.VELOCITY)) >= obj.saccade_detection_threshold),obj.MOV_TYPE ) = obj.SACCADE_TYPE;
% Now we mark every invalid point as noise
            obj.eye_records( (obj.eye_records(:,obj.VALIDITY) == obj.DATA_INVALID),obj.MOV_TYPE ) = obj.NOISE_TYPE;
% This is special case
            obj.eye_records( 1 , obj.MOV_TYPE ) = obj.eye_records( 2 , obj.MOV_TYPE );
% Now we have to expand saccade mark to previous point
            tmp_type = obj.eye_records(2:length(obj.eye_records),obj.MOV_TYPE);
            tmp_type(length(tmp_type)+1) = NaN;
            obj.eye_records( ((obj.eye_records(:,obj.VALIDITY) == obj.DATA_VALID) & (tmp_type(:) == obj.SACCADE_TYPE)),obj.MOV_TYPE) = obj.SACCADE_TYPE;

% Now we check every record to determine its class
            noise_cnt =     sum( obj.eye_records(:,obj.MOV_TYPE)==obj.NOISE_TYPE );
            saccade_cnt =   sum( obj.eye_records(:,obj.MOV_TYPE)==obj.SACCADE_TYPE );
            fixation_cnt =  sum( obj.eye_records(:,obj.MOV_TYPE)==obj.FIXATION_TYPE );

% Estimate probabilities for saccade and fixation
            saccade_probability = saccade_cnt / ( length(obj.eye_records) - noise_cnt );
            fixation_probability = fixation_cnt / ( length(obj.eye_records) - noise_cnt );

% Begin Viterbi algorithm
            for k=1:obj.window_count
% Get begin, end, and length of the current window
                [start_record ended_record] = obj.get_window_parameters(k);
                
% Initializing probabilities
% Declaring states - 1st is the saccade, 2nd is the fixation
                state_l = [1 2];
% Setup starting probabilities
                start_p = [saccade_probability fixation_probability];
% Setup transition probabilities between states - like saccade->fixation or
% fixation->fixation. We assuming that probability to change state is 0.5
                trans_p = [0.95 0.05 ; 0.05 0.95];
% Emissions from states. Watch out! It's a magic numbers!
                emit_p = [0.7 0.2 0 0.1 ; 0.3 0.6 0 0.1];

% Preparing data for Baum-Welch re-estimation algorithm
% This is our pre-estimated observations
                argmax = obj.eye_records( start_record:ended_record,obj.MOV_TYPE );
% And this is probability at the end
                ended_p = [0.1 0.1];
                [start_p, trans_p, emit_p] = obj.baum_welch_reestimation(state_l, argmax, start_p, trans_p, emit_p, ended_p, 0);

% Example from wikipedia and David Conger's code
                T = cell(2);
                for i=1:2
                    T{i} = {start_p(i); state_l(i); start_p(i)};
                end

                for i=start_record:ended_record
                    U=cell(2);
                    for j=1:2
                        total = 0;
                        argmax = [];
                        valmax = 0;
                        for k=1:2
                            Ti = T{k};
                            prob = Ti{1}; v_path = Ti{2}; v_prob=Ti{3};
                            p = emit_p( k,obj.eye_records(i, obj.MOV_TYPE) ) * trans_p(k,j);
                            prob = prob * p;
                            v_prob = v_prob * p;
                            total = total + prob;
                            if( v_prob > valmax )
                                argmax = [v_path; state_l(j)];
                                valmax = v_prob;
                            end
                        end
                        U{j} = {total;argmax;valmax};
                    end
                    T = U;
                end
            
                total = 0;
                argmax = [];
                valmax = 0;
                for i=1:2
                    Ti = T{i};
                    prob = Ti{i}; v_path = Ti{2}; v_prob = Ti{3};
                    total = total + prob;
                    if( v_prob > valmax )
                        argmax = v_path;
                        valmax = v_prob;
                    end
                end
            
% Classify all records based on argmax array (argmax = 2 if saccade and 1
% if fixation)
                L=ended_record - start_record + 1;
                obj.eye_records( start_record:ended_record, obj.MOV_TYPE) =...
                    (     argmax(1:L) - ones( L,1 ) ) * obj.SACCADE_TYPE +...
                    (abs( argmax(1:L) - 2 * ones( L,1 )) )*obj.FIXATION_TYPE;
            end

            obj.eye_records( (obj.eye_records(:,obj.VALIDITY) == obj.DATA_INVALID),obj.MOV_TYPE ) = obj.NOISE_TYPE;
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Complete data classification with I-HMM classifier in :',datestr(now),'\n'));
            end
        end

% Implementation of Baum-Welch algorithm for re-estimating probabilities
        function [start_probability, transition_probability, emission_probability] = ...
                baum_welch_reestimation(obj, ...
                                        states, ...
                                        sample_records, ...
                                        start_probability, ...
                                        transition_probability, ...
                                        emission_probability, ...
                                        ended_probability, ...
                                        recursive_count)            
% Calculate alpha values.
            alpha = zeros(length(sample_records),2);
            T = zeros(2);
            temp_alpha = zeros(2,2);
% For first record in sample records
            for i=1:2
                T(i)=emission_probability(i,sample_records(1)) * start_probability(i);
                alpha(1,i) = T(i);
            end
% And for all others
            for i=2:length(sample_records)
                for j=1:2
                    probability = T(j);
                    for k=1:2
                        p = emission_probability(k, sample_records(i)) * transition_probability(k,j);
                        temp_alpha(k,j) = probability * p;
                    end
                end

                for j=1:2
                    total=0;
                    for k=1:2
                        total=total+temp_alpha(j,k);
                    end
                    T(j)=total;
                    alpha(i,j) = T(j);
                end
            end

% Calculate beta values
            U = zeros(2);
            beta = zeros(length(sample_records),2);
% For first record in sample records
            for i=1:2
                U(i) = ended_probability(i);
                beta(length(sample_records),i) = U(i);
            end
% And for all others
            for i=length(sample_records)-1:-1:1
                for j=1:2
                    total = 0;
                    for k=1:2
                        probability = U(k);
                        p = emission_probability( k, sample_records(i+1)) * transition_probability(k,j) * probability;
                        total = total + p;
                    end
                    if j==1
                        temp1 = total;
                    else
                        temp2 = total;
                    end
                end
                for j=1:2
                    if j==1
                        U(j) = temp1;
                    else
                        U(j) = temp2;
                    end
                    beta (i,j) = U(j);
                end
            end

% Calculate alpha x beta for each state
            alpha_x_beta = zeros(length(sample_records),2);
            for i=1:length(sample_records)
                for j=1:2
                    alpha_x_beta(i,j) = alpha(i,j)* beta(i,j);
                end
            end

% Calculate (alpa(s1) x beta(s1)) + (alpa(s2) x beta(s2)) for each state
            alpha_p_beta = zeros(length(sample_records),1);
            for i=1:length(sample_records)
                total =0;
                for j=1:2
                    total = total + alpha_x_beta(i,j);
                end
                alpha_p_beta(i) = total;
            end
            
            states_probability_total = zeros(2,1);
            states_probability = zeros(length(sample_records),2);
% Calculate probability for each state
            for i=1: length(sample_records)
                for j=1:2
                    states_probability(i,j) = alpha_x_beta(i,j)/alpha_p_beta(i);
                    states_probability_total(j) = states_probability_total(j) + states_probability(i,j);
                end
            end

            probability_given_observation = zeros(length(sample_records),2,3);
            probability_given_observation_total = zeros(2,3);
% Probability that with all observation we can reach each state
            for i=1:length(sample_records)
                for j=1:2
                    switch sample_records(i)
                        case 1  % saccade case
                            probability_given_observation(i,j,1) = states_probability(i,j);
                        case 2  % fixation case
                            probability_given_observation(i,j,2) = states_probability(i,j);
                        case 3  % noise case
                            probability_given_observation(i,j,3) = states_probability(i,j);
                    end
% Calculate total of each column of probabilities
                    for k=1:3
                        probability_given_observation_total(j,k) = probability_given_observation_total(j,k) + ...
                                                                    probability_given_observation(i,j,k);
                    end
                end
            end

% Probability of states transistion given all observations starting from
% 2nd record
            probability_transition_ss = zeros(length(sample_records),1);
            probability_transition_sf = zeros(length(sample_records),1);
            probability_transition_fs = zeros(length(sample_records),1);
            probability_transition_ff = zeros(length(sample_records),1);
            for i=2:length(sample_records)
                probability_transition_ss(i) = (    alpha((i-1),1) * ...
                                                    transition_probability(1,1) * ...
                                                    beta(i,1) * ...
                                                    emission_probability(1,sample_records(i) )) / alpha_p_beta(i);
                probability_transition_sf(i) = (    alpha((i-1),1) * ...
                                                    transition_probability(1,2) * ...
                                                    beta(i,2) * ...
                                                    emission_probability(2, sample_records(i) )) / alpha_p_beta(i);
                probability_transition_fs(i) = (    alpha((i-1),2) * ...
                                                    transition_probability(2,1) * ...
                                                    beta(i,1) * ...
                                                    emission_probability(1,sample_records(i) )) / alpha_p_beta(i);
                probability_transition_ff(i) = (    alpha((i-1),2) * ...
                                                    transition_probability(2,2) * ...
                                                    beta(i,2) * ...
                                                    emission_probability(2,sample_records(i) )) / alpha_p_beta(i);
            end
            probability_states_transition_ss_total = sum(probability_transition_ss);
            probability_states_transition_sf_total = sum(probability_transition_sf);
            probability_states_transition_fs_total = sum(probability_transition_fs);
            probability_states_transition_ff_total = sum(probability_transition_ff);

% Emission probability re-estimation
            for i=1:2
                for k=1:3
                    emission_probability(i,k) = probability_given_observation_total(i,k) / states_probability_total(i);
                end
            end
            
% Transition probability re-estimation
            transition_probability(1,1) = probability_states_transition_ss_total / states_probability_total(1);
            transition_probability(1,2) = probability_states_transition_sf_total / states_probability_total(1);
            transition_probability(2,1) = probability_states_transition_fs_total / states_probability_total(2);
            transition_probability(2,2) = probability_states_transition_ff_total / states_probability_total(2);

% Start probability re-estimation
            for i=1:2
                start_probability(1,i) = states_probability(1,i);
            end

% End probability reestimation
            for i=1:2
                ended_probability(1,i) = states_probability(length(sample_records),i)/states_probability_total(i);
            end
            
            if (recursive_count < obj.baum_welch_reiterations_count)
                [start_probability, transition_probability, emission_probability] = ...
                    obj.baum_welch_reestimation(states, sample_records, start_probability,transition_probability, emission_probability, ended_probability, recursive_count+1);
            end
        end

% Access interface for class properties
        function set.saccade_detection_threshold(obj, value)
            obj.saccade_detection_threshold = value;
        end

        function set.viterbi_sample_size(obj, value)
            obj.viterbi_sample_size = value;
        end

        function set.baum_welch_reiterations_count(obj, value)
            obj.baum_welch_reiterations_count = value;
        end

        function result = get.window_count(obj)
            result = ceil(length (obj.eye_records) / obj.viterbi_sample_size);
        end

        function [start_window end_window] = get_window_parameters(obj, count)
            start_window = obj.viterbi_sample_size * (count-1) + 1;
            end_window = start_window + obj.viterbi_sample_size - 1;
            if (end_window > length(obj.eye_records) )
                end_window = length( obj.eye_records);
            end
        end

    end

end
