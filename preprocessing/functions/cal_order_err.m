function edf = cal_order_err(edf)
% Calculate angular and order error
            % retrieve target position in dva
            tarx = edf.param.tarx_deg(1:16,:);
            tary = edf.param.tary_deg(1:16,:);
            
            % retrieve primary saccade and the closest saccade endpoint
            % position
            psacx = edf.cal.primary_sac_xend;
            psacy = edf.cal.primary_sac_yend;
            asacx = edf.cal.acc_sac_tar_xend;
            asacy = edf.cal.acc_sac_tar_yend;
            % psactx = data1{ii,jj}.cal.primary_sac_tar_xend; % primary
            % saccade for the target
            % psacty = data1{ii,jj}.cal.primary_sac_tar_yend;
            
            % calculate the polar angles for the target
            [tar_theta,~] = cart2pol(tarx,tary); tar_theta = rad2deg(tar_theta);
            % polar angle for the primary saccade (to any stimulus)
            [psac_theta,~] = cart2pol(psacx,psacy); psac_theta = rad2deg(psac_theta);
            % polar angle for the closest saccade to the target
            [asac_theta,~] = cart2pol(asacx,asacy); asac_theta = rad2deg(asac_theta);
            %[psact_theta,psact_rho] = cart2pol(psactx,psacty); psact_theta = rad2deg(psact_theta);
              
            % calculate polar angle differences between saccade and target
            % (angular error)
            % primary saccade minus target angle
            p_minus_tar_theta = psac_theta -tar_theta;
            % the closest saccade minus target angle
            a_minus_tar_theta = asac_theta -tar_theta;
            % pt_minus_tar_theta = psact_theta -tar_theta;
            
            % if the angular different is more than 180 degrees
            % that means the saccade is more than 180 degrees clockwise
            % relative to the target
            % that is, the saccade is less than 180 degrees counterclockwise
            % relative to the target -> thus we will substract 360 from the
            % original angular error (200 degrees -> -160 degrees)
            ind1 = find(p_minus_tar_theta>180);
            p_minus_tar_theta(ind1) = p_minus_tar_theta(ind1)-360;
            ind2 = find(a_minus_tar_theta>180);
            a_minus_tar_theta(ind2) = a_minus_tar_theta(ind2)-360;

            % similarly, if the angular difference is less than -180
            % degrees
            % that means the saccade is more than 180 degrees
            % counterclockwise to the target
            % that is, the saccade is less than or EQUAL TO 180 degrees clockwise to
            % the target -> thus we will add 360 to the original angular
            % error (-200 degrees -> 160 degrees)
            % note that we have <= -180 here, meaning that if the angular
            % error is -180, it will be converted to 180. There will only
            % be 180 degrees, not -180 degrees, in our calculation
            clear ind1 ind2
            ind1 = find(p_minus_tar_theta<=-180);
            p_minus_tar_theta(ind1) = p_minus_tar_theta(ind1)+360;
            ind2 = find(a_minus_tar_theta<=-180);
            a_minus_tar_theta(ind2) = a_minus_tar_theta(ind2)+360;
            
            %             % plot and check - this is for debugging
            %             figure(1);clf
            %             subplot(1,2,1)
            %             polarplot(tar_theta,tar_rho,'*')
            %             hold on
            %             subplot(1,2,2)
            %             scatter(tarx,tary,'*')
            
            % Let's further calculate target and saccade quadrant based on
            % its polar angle
            % 0 <= error < 90 degrees: first quadrant
            % 90 <= error < 180 degrees: seccond quadrant
            % error == 180 or -180 < error < -90: third quadrant
            % -90 <= error < 0: fourth quadrant
            % Note that we only have 180, not -180 degrees using our
            % calculation above

            % target quadrant
            tar_quad = tar_theta;
            tar_quad(tar_theta>=0 & tar_theta <90) = 1;
            tar_quad(tar_theta>=90 & tar_theta <180) = 2;
            tar_quad(tar_theta == 180 | (tar_theta>-180 & tar_theta <-90)) = 3;
            tar_quad(tar_theta>=-90 & tar_theta <0) = 4;
            
            % primary saccade quadrant
            psac_quad = psac_theta;
            psac_quad(psac_theta>=0 & psac_theta <90) = 1;
            psac_quad(psac_theta>=90 & psac_theta <180) = 2;
            psac_quad(psac_theta == 180 | (psac_theta>-180 & psac_theta <-90)) = 3;
            psac_quad(psac_theta>=-90 & psac_theta <0) = 4;
            
            % the closest saccade to the target
            asac_quad = asac_theta;
            asac_quad(asac_theta>=0 & asac_theta <90) = 1;
            asac_quad(asac_theta>=90 & asac_theta <180) = 2;
            asac_quad(asac_theta == 180 |(asac_theta>-180 & asac_theta <-90)) = 3;
            asac_quad(asac_theta>=-90 & asac_theta <0) = 4;
            
            % calculate quadrant difference between saccade and target
            p_minus_tar_quad = psac_quad - tar_quad; % primary saccade minus target
            a_minus_tar_quad = asac_quad - tar_quad; % closest saccade minus target
            
            % store calculated angular, and quadrant errors
            edf.cal.psac_theta = psac_theta;
            edf.cal.asac_theta = asac_theta;
            edf.cal.tar_theta = tar_theta;
            edf.cal.p_minus_tar_theta = p_minus_tar_theta;
            edf.cal.a_minus_tar_theta = a_minus_tar_theta;
            
            edf.cal.psac_quad = psac_quad;
            edf.cal.asac_quad = asac_quad;
            edf.cal.tar_quad = tar_quad;
            edf.cal.p_minus_tar_quad = p_minus_tar_quad;
            edf.cal.a_minus_tar_quad = a_minus_tar_quad;      
            
            % Now let's calculate order transposition error
            % retrieve the orders
            % order of the target (which item was probed during that trial)
            tar_ord = edf.param.probe_order(1:16,:);
            % the order that each quadrant appeared
            % each row is one trial
            % each column represent the ordinal rank, value represents the quadrant that appears for that rank
            quad_ord = edf.param.quad_order(1:16,:); 

            % let's construct a new variable, ord_quad, such that  each 
            % column represent the quadarnt number, value represents the 
            % ordinal rank that quadrant appears
            for kk = 1:16 % for each trial
                ord_quad(kk,1) = find(quad_ord(kk,:)==1);
                ord_quad(kk,2) = find(quad_ord(kk,:)==2);
                ord_quad(kk,3) = find(quad_ord(kk,:)==3);
                ord_quad(kk,4) = find(quad_ord(kk,:)==4);
            end
            
            % which item did the saccade respond to?
            for kk = 1:16 % for each trial
                if ~isnan(psac_quad(kk)) % if there is a primary saccade detected for that trial
                    psac_ord(kk) = ord_quad(kk,psac_quad(kk)); % let's calculate the serial position (which item) that the saccade responded to   
                else
                    psac_ord(kk) = nan;
                    
                end   
                % similary for the closest saccade to the target
                if ~isnan(asac_quad(kk))
                    asac_ord(kk) = ord_quad(kk,asac_quad(kk));
                else
                    asac_ord(kk) = nan;
                end
            end

            % transposition error is calculated as the serial position that
            % the saccade responded to, minus the serial position of the
            % actual target
            % if this error is +1, it means that this saccade responds to
            % one item later to the target 
            % if this error is -2, it means that this saccade responds to
            % two items earlier to the target
            p_minus_tar_ord = psac_ord' - tar_ord;
            a_minus_tar_ord = asac_ord' - tar_ord;
            
            % store the transposition error and serial positions
            edf.cal.psac_ord = psac_ord'; % which item primary saccade responded to
            edf.cal.asac_ord = asac_ord'; % which item closest saccade responded to
            edf.cal.tar_ord = tar_ord; % which item is the target
            edf.cal.p_minus_tar_ord = p_minus_tar_ord; % transposition error for primary saccade
            edf.cal.a_minus_tar_ord = a_minus_tar_ord; % transposition error for the closest saccade


end

