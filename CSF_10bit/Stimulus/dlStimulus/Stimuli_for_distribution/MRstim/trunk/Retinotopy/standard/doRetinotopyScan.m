function doRetinotopyScan(params)
% doRetinotopyScan - runs retinotopy scans
%
% doRetinotopyScan(params)
%
% Runs any of several retinotopy scans
%
% 99.08.12 RFD wrote it, consolidating several variants of retinotopy scan code.
% 05.06.09 SOD modified for OSX, lots of changes.

% defaults
if ~exist('params', 'var')
	error('No parameters specified!');
end

% quit key
try, 
    quitProgKey = params.display.quitProgKey;
catch,
    quitProgKey = KbName('q');
end;

% make/load stimulus
switch params.experiment
    case {'Translating Bars 8 Pass (Dumoulin)','Translating Bars 8 Pass (contrast)'},
        [stimulus]= makeRetinotopyStimulus_bars8Pass(params);
    otherwise,
        stimulus = makeRetinotopyStimulus(params);
end;



% allow stimulus to start a little earlier
% stimulus.offset = 3; % secs
% if stimulus.offset>0
%     stimulus.unclipped = stimulus; % keep the old one just in case
%     
%     stimulus.seqtimging = stimulus.seqtiming-stimulus.offset;
%     keep = stimulus.seqtiming>=0;
%     stimulus.seq = stimulus.seq(keep);
%     stimulus.seqtiming = stimulus.seqtiming(keep);
%     stimulus.fixSeq = stimulus.fixSeq(keep);
%     if exist('otherSequence')
%         otherSequence = otherSequence(keep);
%     end
% end

% loading mex functions for the first time can be
% extremely slow (seconds!), so we want to make sure that 
% the ones we are using are loaded.
KbCheck;GetSecs;WaitSecs(0.001);

try,
    % check for OpenGL
    AssertOpenGL;
    
    % to skip annoying warning message on display (but not terminal)
    Screen('Preference','SkipSyncTests', 1);
    
    % Open the screen
    params.display                = openScreen(params.display);
    params.display.devices        = params.devices;

    % to allow blending
    Screen('BlendFunction', params.display.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Store the images in textures

    stimulus = createTextures(params.display,stimulus,0);

    % getting to the source rather than doScan->doTrial->showStimulus 
    for n = 1:params.repetitions,
        % set priority
        Priority(params.runPriority);
        
        % reset colormap?
        resetColorMap(params);
        
        % flush the device
        deviceUMC('read',params.devices.UMCport);
        
        % wait for go signal
        dispStringInCenter(params.display,'Waiting for scanner to start',0.6,[],'white');
%         tmp = Screen('MakeTexture',params.display.windowPtr, ...
%             zeros(768,1024));
%         Screen('DrawTexture', params.display.windowPtr, tmp);
%         Screen('Flip',params.display.windowPtr);

        [junk time0]=deviceUMC('wait for trigger',params.devices.UMCport);  
        %disp(sprintf('[%s]:%s',mfilename,'Waiting for scanner to start'));

        % go
        [response, timing, quitProg] = showScanStimulus(params.display,stimulus,time0);

        
        % reset priority
        Priority(0);
        
        % get performance
        [pc,rc] = getFixationPerformance(params.fix,stimulus,response);
        disp(sprintf('[%s]:Fixation dot task: percent correct: %.1f %%, reaction time: %.1f secs',mfilename,pc,rc));
        
        % get other detection performance
        if exist('otherSequence') && ~isempty(otherSequence),
            [pc,rc,nn] = getDetectionPerformance(params.fix,stimulus,response,otherSequence);
             disp(sprintf('[%s]:Other task: percent correct: %.1f %%, reaction time: %.1f secs',mfilename,pc,rc));
        end;

        
        % save 
        if params.savestimparams,
            filename = ['~/Desktop/' datestr(now,30) '.mat'];
            save(filename);                % save parameters
            disp(sprintf('[%s]:Saving in %s.',mfilename,filename));
        end;
        
        % keep going?
        if quitProg, % don't keep going if quit signal is given
            break;
        end;
    end;
    
    % Close the one on-screen and many off-screen windows
    closeScreen(params.display);
catch,
    % clean up if error occurred
    Screen('CloseAll');
    setGamma(0);
    Priority(0);
    ShowCursor;
    rethrow(lasterror);
end;


return;


function resetColorMap(params);
switch params.experiment,
    case {'8 bars (slow)'},
        % get subject input on stimulus contrast
        answer = 'a';
        while (~isnumeric(answer) || isempty(answer)),
            % say it just in case the screen is in mirror mode and you
            % cannot see the screen
            if params.display.screenNumber == 0,
                eval('system(sprintf(''say please enter percent stimulus contrast.''));'); 
            end;
            answer = input('Please enter percent stimulus contrast [100]: ');
        end;
        sz = size(params.display.gamma,1)*([-1 1]*answer/100/2+0.5);
        if sz(1)==0,sz(1)=1;end;
        % make new gamma
        putgamma = zeros(256,3);
        putgamma(2:255,:) = params.display.gamma(round(linspace(sz(1),sz(2),254)),:);
        % load gamma
        putgamma(1,:)   = params.display.fixColorRgb(1,1:3)./255;
        putgamma(256,:) = params.display.fixColorRgb(2,1:3)./255;
        Screen('LoadNormalizedGammaTable', params.display.screenNumber,putgamma);
        
    case {'8 bars (LMS)','8 bars (LMS) with blanks'},
        % get subject input on stimulus type 
        answer = 'a';
        while (~isnumeric(answer) || isempty(answer) || answer>3 || answer<1),
            % say it just in case the screen is in mirror mode and you
            % cannot see the screen
            if params.display.screenNumber == 0,
                eval('system(sprintf(''say please enter  stimulus type 1, 2, or 3.''));'); 
            end;
            answer = input('Please enter  stimulus type [1=LMS,2=L-M,3=S]: ');
        end;
        if answer == 1,
            stimtype = [ 1 1 1];
        elseif answer == 2,
            stimtype = [-1 1 0];
        else,
            stimtype = [0 0 1];
        end;
        
        % get subject input on stimulus contrast
        answer = 'a';
        while (~isnumeric(answer) || isempty(answer) || answer<0 || answer>100),
            % say it just in case the screen is in mirror mode and you
            % cannot see the screen
            if params.display.screenNumber == 0,
                eval('system(sprintf(''say please enter percent stimulus contrast 0 to 100.''));'); 
            end;
            answer = input('Please enter percent stimulus contrast: ');
        end;
        stimcontrast = answer;
        
        
        sz = size(params.display.gamma,1)*([-1 1]*answer/100/2+0.5);
        if sz(1)==0,sz(1)=1;end;
        % make new gamma        
        newgamma = create_LMScmap(params.display,stimtype.*(stimcontrast./100));
        putgamma = zeros(256,3);
        if size(newgamma,1)~=256,
            putgamma(2:255,:) = newgamma(round(linspace(1,size(newgamma,1),254)),:);
        end;
        % load gamma
        putgamma(1,:)   = params.display.fixColorRgb(1,1:3)./255;
        putgamma(256,:) = params.display.fixColorRgb(2,1:3)./255;
        Screen('LoadNormalizedGammaTable', params.display.screenNumber,putgamma);

        
        
    otherwise,
        % well.. nothing
end;
return;
