function edf = detect_epoch_exp1_exp2(edf,set)
% Detect the onset and offset of a trial-related message
% Including the start and end of a trial
% And the start and end of a task epoch

msg = set.msg; % prototype messages in a trial

% index of the messages that contains the task epoch message info
ind = find(contains(edf.default_events.Messages.info,msg));
% content of the messages
txt = {edf.default_events.Messages.info{ind}};
% convert text to numbers for messages for easier analysis in the future
% the first to the last message in set.msg corresponds to a number from 1
% to length(set.msg)
txt_num = zeros(size(ind));
for ii = 1:length(set.msg) % for each message
    txt_num(contains(txt,set.msg{ii})) = ii;
end
% timing of the messages
tm = edf.default_events.Messages.time(ind);

% Let's find out all the trial IDs
srt_idx = find(txt_num == 1);
id = split(txt(srt_idx),' '); id = id(:,:,2); id = cellfun(@(c) str2num(c),id); % trial IDs
id = id + 1; % we want to add 1 to the ID because the ID starts from 0
end_idx = find(txt_num == length(set.msg));
% remove the 17th trial as it is a baseline trial
if find(id == 17) % if we found a 17th trial
    id(id==17) = []; srt_idx(id==17) == []; end_idx(id==17) == []; % remove that trial from the id and indexes
end

msgText = nan(16,length(set.msg)); % actual content of the message
msgTime = nan(16,length(set.msg)); % timing of the message
% for each trial ID, let's find out the corresponding messages between the
% current trial starts and current trial end
for ii = 1:length(id)
    curr_srt = srt_idx(ii); % message index of the current trial start (relative to txt)
    curr_end = end_idx(ii); % message index of the current trial end (relative to txt)
    
    % What messages do we have in between trial start and trial end?
    
    if curr_srt == curr_end - 1 | tm(curr_end)-tm(curr_srt) > 25000 | tm(curr_end)-tm(curr_srt) < 10000 
        % if there is no messages in between
        % of if current trial length is longer than 25 s or shorter than
        % 10 s
        % we will skip this trial
        % because there is no way to calculate the event timing based on trial start and end
        % or because the trial is too long or too short, so there must be
        % something wrong with this trial
        
        % we still want to store the start and the end of that trial
         msgText(id(ii),[1,length(set.msg)]) = txt_num([curr_srt,curr_end]);
         msgTime(id(ii),[1,length(set.msg)]) = tm([curr_srt,curr_end]);

        continue;
        warning(['Skip trial ' num2str(id(ii))])
        
    else % if there is at least one message in between
        % calculate what messages are there
        msg_have = txt_num(curr_srt+1 : curr_end-1);
        % what messages should be there
        msg_should = 2:length(set.msg)-1;
        
        % if we have all the messages recorded
        if isequal(msg_have,msg_should)
            % store the message content and time
            msgText(id(ii),:) = txt_num(curr_srt:curr_end); % the content is converted to numbers 1-length(set.msg), corresponding to each msg
            msgTime(id(ii),:) = tm(curr_srt:curr_end); % message time (relative to sample time)
            
        else % we miss some messsages in between, but still have at least one message recorded
            warning(['Miss some messages for trial ' num2str(id(ii))])
            
            % let's record the timing of the messages that we have first
            msgText(id(ii),[1,msg_have,length(set.msg)]) = txt_num([curr_srt:curr_end]);
            msgTime(id(ii),[1,msg_have,length(set.msg)]) = tm([curr_srt:curr_end]);
            
            % Then we will calculate time of the missed messages based on the
            % message that we have and the planned timing of the events
            
            % What is the planned time of all the events?
            planned_tm = edf.param.time; % output from load_params
            % from left to right should be the END time of the 2nd to the second to last
            % message (event), in a cumulative fashion
            % Note that in our current task, we have iti at the beginning of
            % each trial, so the last column have a smaller time than all the
            % rest of the events
            % e.g., 3500, 4000, ...., 2000
            
            % we want to transform this timing matrix such that it
            % marks the START, not the end time of each event
            % To do that, we shift the columns to the right by 1
            iti_end_tm = planned_tm(:,end);
            planned_tm(:,2:end) = planned_tm(:,1:end-1);
            planned_tm(:,1) = iti_end_tm;

            % After shifting the columns, each cell represents the ONSET of 
            % a specific event right now, from fixation to iti
            
            % to make it consistent across the msgTime, we will add one column
            % at the beginning and at the end of this matrix respectively, to
            % account for timing of the first (trial start) and the last (trial
            % end) event
            planned_tm = [zeros(size(edf.param.time,1),1) planned_tm zeros(size(edf.param.time,1),1)];
            
            % What messages do we miss?
            msg_miss = setdiff(msg_should,msg_have);
            
            clear msg_miss_tm msg_close close_idx msg_close_tm miss_minus_close
            % Let's calculate the timing of the missed messages based on the
            % closest message that we have
            for mm = 1:length(msg_miss) % for each missed message
                % what is the closest message that we have?
                [~,close_idx] = min(abs(msg_miss(mm)-msg_have)); % close_idx, the location of the item in msg_have that is closest to the missed messsage
                msg_close = msg_have(close_idx); % closest message
                
                % what is the time interval between the closest message that we
                % have and the missed message?
                % planned missed message time minus planned recorded message
                % time
                msg_close_tm = tm(curr_srt+close_idx);
                miss_minus_close = planned_tm(id(ii),msg_miss(mm)) - planned_tm(id(ii),msg_close);
                
                % calculate missed message time
                msg_miss_tm(mm) = msg_close_tm + miss_minus_close;
            end
            
            % store the missed message text and time
            msgText(id(ii),msg_miss) = msg_miss;
            msgTime(id(ii),msg_miss) = msg_miss_tm;
        end
    end
    
end

    % We want to retrieve the index of the sample when each message appears
    [~,msg_srt_idx] = arrayfun(@(t) min(abs(t-edf.samples.time)),msgTime);
    msg_srt_idx(isnan(msgTime)) = nan;
    
    % to get the index of the end of the message, we just need to shift the
    % columns to the left by 1
    % and remember for the last column (end of the trial), its end is the
    % start of the next trial
    msg_end = msg_srt_idx(:,1);
    msg_end_idx(:,1:length(set.msg)-1) = msg_srt_idx(:,2:end);
    msg_end_idx(:,length(set.msg)) = [msg_end(2:end);length(edf.samples.time)];
    
    edf.events.msg.txt = msgText;
    edf.events.msg.time = msgTime;
    edf.events.msg.ind_srt = msg_srt_idx;
    edf.events.msg.ind_end = msg_end_idx;
      
    % store trial and epoch info in edf.samples
    edf.samples.trial = zeros(size(edf.samples.time));
    edf.samples.msg = edf.samples.trial;
    for ii = 1:size(msgText,1) % for each trial
    edf.samples.trial(msg_srt_idx(ii,1):msg_end_idx(ii,length(set.msg))) = ii;
    for jj = 1:size(msgText,2) % for each event
        if ~isnan(msg_srt_idx(ii,jj)) & ~isnan(msg_end_idx(ii,jj))
    edf.samples.msg(msg_srt_idx(ii,jj):msg_end_idx(ii,jj)) = msgText(ii,jj);
        end
    end
    end

% We double check the time
% If the timing differences between events

% d = diff(msgTime,1,2);

% 
% % Let's store the onset and offset of each trial and epoch separately
% edf.events.num_trial = ntr;
% if mod(length(msgText),ntr) == 0
% edf.events.msg.txt= reshape(msgText,ntr,[]);
% edf.events.msg.time = reshape(msgTime,ntr,[]);
% else % cannot be divided
%     for tt = 1:length(set.msg)
%         ind_txt = find(msgText==tt);
%         if length(ind_txt) == ntr
%         edf.events.msg.txt(:,tt) = msgText(ind_txt);
%         edf.events.msg.time(:,tt) = msgTime(ind_txt);
%         elseif length(ind_txt) > ntr
%         edf.events.msg.txt(:,tt) = msgText(ind_txt(1:ntr));
%         edf.events.msg.time(:,tt) = msgTime(ind_txt(1:ntr));
% 
%         else
%         edf.events.msg.txt(:,tt) = [msgText(ind_txt);nan(ntr-length(ind_txt),1)];
%         edf.events.msg.time(:,tt) = [msgTime(ind_txt);nan(ntr-length(ind_txt),1)];
%         end
%     end
% end
% 
% % onset and offset of a trial
% msgTime = edf.events.msg.time;
% msgText = edf.events.msg.txt;
% for ii = 1:ntr
% for jj = 1:size(msgText,2)
% [val,ind] = min(abs(edf.samples.time-msgTime(ii,jj)));
% edf.events.msg.ind_srt(ii,jj) = ind(1);
% if size(msgText,2) > 1 % more than 1 message
% if jj == 1 % the first message (e.g., start of a trial)
%     if ii == 1 % the first trial
%         edf.events.msg.ind_end(ii,jj) = nan;
%     else % not the first trial
%     edf.events.msg.ind_end(ii-1,size(msgText,2)) = ind(1); % end of the last message, last trial
%     end
% else % not the first message
%     if ii == ntr & jj == size(msgText,2) % the last message and the last trial
% edf.events.msg.ind_end(ii,jj) = length(edf.samples.time); % last time point
%     else
% edf.events.msg.ind_end(ii,jj-1) = ind(1); % end of the last epoch in the same trial
%     end
% end
% elseif size(msgText,2) == 1 % only has 1 message
%         if ii == 1 % the first trial
%         edf.events.msg.ind_end(ii,jj) = nan;
%         elseif ii == ntr % the last trial
%             edf.events.msg.ind_end(ii,jj) = length(edf.samples.time); % last time point
%         else % in between
%            edf.events.msg.ind_end(ii-1,jj) = ind(1); 
%         end
% end
% 
% end
% end


end

