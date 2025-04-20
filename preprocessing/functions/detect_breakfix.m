function edf = detect_breakfix(edf)
% Detect trials that break fixation
            % endpoint fixation position of saccade
            end_pos = sqrt(edf.events.sac_dc.endfix_avg_x.^2 + edf.events.sac_dc.endfix_avg_y.^2)';
            % messages during the saccade offset
            end_msg = edf.events.sac_dc.msg_end;
            % trials to which the saccade belongs
            end_trial = edf.events.sac_dc.trial;
            
            % initiate a field called 'breakfix' under 'events.sac_dc' to
            % store saccade indexes where subjects broke fixation (move their
            % eyes more than 2 degrees from the center of the screen)
            edf.events.sac_dc.breakfix = zeros(size(end_trial));
            
            % another field 'breakfix_trial' under edf.trackloss to store the corresponding
            % trial number
            breakfix_trial = [];
            
            for tt = 1:edf.samples.ntrial % 16 trials
                % find row indexes corresponding to a specific trial,
                % with messages 2-9, and with its saccade end position more
                % than 3 degrees from the screen center
                rows = find(end_trial==tt & ismember(end_msg,[3:10]) & end_pos > 2);
                if ~isempty(rows)
                    edf.events.sac_dc.breakfix(rows) = 1;
                    breakfix_trial = [breakfix_trial;tt];
                end
            end    
 
            edf.trackloss.breakfix = breakfix_trial;
            edf.trackloss.nanerr = find(isnan(edf.cal.primary_sac_err)); % if there is no saccade that enters the AOI of ANY stimulus
            edf.trackloss.rmtrial = unique([edf.trackloss.breakfix;edf.trackloss.nanerr;edf.trackloss.nodc_trial']);

end

