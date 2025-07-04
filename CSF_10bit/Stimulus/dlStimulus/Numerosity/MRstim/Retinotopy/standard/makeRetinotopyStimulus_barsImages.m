function [stimulus] = makeRetinotopyStimulus_barsImages(params)
% makeRetinotopyStimulus - make various retinotopy stimuli
%
% stimulus = makeRetinotopyStimulus_bars(params)
%
% Matlab code to generate various retinotopy stimuli
% Generates one full cycle, as well as the sequence for the entire scan.
%
% 99.09.15 RFD: I fixed the sequence generation algorithm so that
%   timing is now frame-accurate.  The algorithm now keeps track
%   of timing error that accumulates due to rounding to the nearest
%   frame and corrects for that error when it gets to be more than 
%   half a frame.  
%   The algorithm also randomely reverses the drift direction, rather
%   than reversing every half-an image duration.
% 2005.06.15 SOD: changed for OSX - stimulus presentation will now be 
%                 time-based rather than frame based. Because of bugs
%                 with framerate estimations.

% 2009.09.27 BMH: Set up for contour integration stimuli

%Contour stimulus parameters
%contourOrientation=params.contour.contourOrientation; %Contour orientation relative to bars. 0 is parrallel to bars
%contourSF=1/3; %Contour spatial frequence, in degrees (1/3=3 cycles per degree)
%contourBandpass=params.contour.contourBandpass; %Contour orientation bandpass, in degrees
pauseDuration=4;
softEdgeWindow=0; %

flipUpDown=0;

% load('/Users/student/Documents/Martijn/retintopimg.mat','-mat');
% 
% disp(sprintf('[%s]:resizing images to 768x768.',mfilename));
% for i = 1:length(naturalimg)
%     naturalimg(i).image = imresize(mat2gray(naturalimg(i).image,[0
%     255]),[768 768]);
% end
% 
% save('/Users/student/Documents/Martijn/retintopimg(768).mat','naturalimg');

global retinotopimg;

if isempty(retinotopimg)
    disp(sprintf('[%s]: loading natural images (538x538).',mfilename));
    load('/Users/lab/Documents/MATLAB/MRstim/trunk/Retinotopy/retinotopimg(538).mat','-mat');
end


% various time measurements:
duration.stimframe          = 1./params.temporal.frequency./params.temporal.motionSteps;
duration.scan.seconds       = params.ncycles*params.period;
duration.scan.stimframes    = params.ncycles*params.period./duration.stimframe;
duration.cycle.seconds      = params.period;
duration.cycle.stimframes   = params.period./duration.stimframe;
duration.prescan.seconds    = params.prescanDuration;
duration.prescan.stimframes = params.prescanDuration./duration.stimframe;


% load matrix or make it
if ~isempty(params.loadMatrix),
    % we should really put some checks that the matrix loaded is
    % appropriate etc.
    load(params.loadMatrix);
    halfNumImages = params.numImages./2;
    disp(sprintf('[%s]:loading images from %s.',mfilename,params.loadMatrix));
%    disp(sprintf('[%s]:size stimulus: %dx%d pixels.',mfilename,n,m));
else
   outerRad = params.radius;
    innerRad = params.innerRad;
    wedgeWidth = params.wedgeWidth;
    ringWidth = params.ringWidth;

    halfNumImages = params.numImages;
    numMotSteps = params.temporal.motionSteps;
    numSubRings = params.numSubRings;
    numSubWedges = params.numSubWedges;

    %%% Set check colormap indices %%%
    %bk = findName(params.display.reservedColor,'background');
    %minCmapVal = max([params.display.reservedColor(:).fbVal])+1;
    %maxCmapVal = params.display.numColors-1;
    bk = params.display.backColorIndex;
    
    
    minCmapVal = min([params.display.stimRgbRange]);
    maxCmapVal = max([params.display.stimRgbRange]);


    %%% Initialize image template %%%
    m=angle2pix(params.display,2*outerRad); 
    n=angle2pix(params.display,2*outerRad);
    
    % should really do something more intelligent, like outerRad-fix
    switch(lower(params.display.fixType))
        case 'left disk',
            [x,y]=meshgrid(linspace( 0,outerRad*2,n),linspace(outerRad,-outerRad,m));
            outerRad = outerRad.*2;
        case 'right disk',
            [x,y]=meshgrid(linspace(-outerRad*2,0,n),linspace(outerRad,-outerRad,m));
            outerRad = outerRad.*2;
        otherwise,
            [x,y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));
    end;
    
    % here we crop the image if it is larger than the screen 
    % seems that you have to have a square matrix, bug either in my or
    % psychtoolbox' code - so we make it square
    if m>params.display.numPixels(2),
        start  = round((m-params.display.numPixels(2))/2);
        len    = params.display.numPixels(2);
        y = y(start+1:start+len, start+1:start+len);
        x = x(start+1:start+len, start+1:start+len);
        m = len;
        n = len;
    end;
    disp(sprintf('[%s]:size stimulus: %dx%d pixels.',mfilename,n,m));

    % r = eccentricity; theta = polar angle
    r = sqrt (x.^2  + y.^2);
    theta = atan2 (y, x);					% atan2 returns values between -pi and pi
    theta(theta<0) = theta(theta<0)+2*pi;	% correct range to be between 0 and 2*pi

    
    % loop over different orientations and make checkerboard
    % first define which orientations
    orientations = (0:45:360)./360*(2*pi); % degrees -> rad
    orientations = orientations([3 2 5 4 7 6 1 8]);  %[6 3 8 5 2 7 4 1]);
    remake_xy    = zeros(1,params.numImages)-1;
    remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
    original_x   = x;
    original_y   = y;
    % step size of the bar
    step_nx      = duration.cycle.seconds./params.tr/8;
    step_x       = (2*outerRad) ./ step_nx;
    step_startx  = (step_nx-1)./2.*-step_x - (ringWidth./2);
    %[0:step_nx-1].*step_x+step_startx+ringWidth./2
    disp(sprintf('[%s]:stepsize: %f degrees.',mfilename,step_x));
    
    % if we create colored bars we want to make the edges soft.
    switch params.experiment,
        case {'8 bars (LMS) with blanks'}
            edgewidth = 25; % .25cpd on a 600 pixel(3deg screen)
            softmask  = makecircle(m-2*edgewidth,m,edgewidth);
        otherwise,
            softmask = ones(m);
    end;

    % Loop that creates the final images
    fprintf('[%s]:Creating %d images:',mfilename,params.numImages);
    images=zeros(m,n,params.numImages*params.temporal.motionSteps,'uint8');
    for imgNum=1:params.numImages
        
        if remake_xy(imgNum) >=0,
            x = original_x .* cos(remake_xy(imgNum)) - original_y .* sin(remake_xy(imgNum));
            y = original_x .* sin(remake_xy(imgNum)) + original_y .* cos(remake_xy(imgNum));

            % Calculate checkerboard.
            % Wedges alternating between -1 and 1 within stimulus window.
            % The computational contortions are to avoid sign=0 for sin zero-crossings

            %syntax for calling synthetic_contours: im=synthetic_contours(sizeInPix,sizeInDeg,contourWidthInDeg, orientationInDeg, orientationBPinDeg);

            % reset starting point
            loX = step_startx - step_x;
  
        end;
        
        prevInd = [0 0 0 0];
        
        for ii=1:3,
            
            rindex = ceil( rand() * length( retinotopimg ) );
            
            while ~isempty( intersect( rindex, prevInd ) )
                rindex = ceil( rand() * length( retinotopimg ) );
            end
            
            prevInd(ii) = rindex;
            if flipUpDown==1
                unscaledContours(:,:,ii) = flipud(retinotopimg(rindex).image);
            else
                unscaledContours(:,:,ii) = retinotopimg(rindex).image;
            end
        end;


        switch params.type;
            case 'bar'
                loEcc = innerRad;
                hiEcc = outerRad;
                loX   = loX + step_x;
                hiX   = loX + ringWidth;
            otherwise,
                error('Unknown stimulus type!');

        end
        % This isn't as bad as it looks
        % Can fiddle with this to clip the edges of an expanding ring - want the ring to completely 
        % disappear from view before it re-appears again in the middle.

        % Can we do this just be removing the second | from the window
        % expression? so...
        if softEdgeWindow>0
            loXin=loX+(softEdgeWindow);
            hiXin=hiX-(softEdgeWindow);
            outerRadIn=outerRad-(softEdgeWindow);
            loXout=loXin-softEdgeWindow;
            hiXout=hiXin+softEdgeWindow;
            outerRadOut=outerRadIn+softEdgeWindow;
            window=( (x>=loXout & x<=hiXout) & r<=outerRadOut);
            tmpvar=zeros(m,n);
            tmpvar(window) = 1;
            loXfade=(x<loXin & x>loXout).*(cos((x-loXin)*pi/softEdgeWindow)./2+.5);
            hiXfade=(x>hiXin & x<hiXout).*(cos((hiXin-x)*pi/softEdgeWindow)./2+.5);
            outerRadFade=(r>outerRadIn & r<outerRadOut).*(cos((r-outerRadIn)*pi/softEdgeWindow)./2+.5);
            tmpvar(loXfade>0)=tmpvar(loXfade>0).*loXfade(loXfade>0);
            tmpvar(hiXfade>0)=tmpvar(hiXfade>0).*hiXfade(hiXfade>0);
            tmpvar(outerRadFade>0)=tmpvar(outerRadFade>0).*outerRadFade(outerRadFade>0);            
            for jj=1:size(unscaledContours,3)
                checks(:,:,jj)=unscaledContours(:,:,jj).*tmpvar;
                img(:,:,jj)=minCmapVal+ceil((maxCmapVal-minCmapVal) * ((checks(:,:,jj)+1)./2));
            end
        else
            checks = unscaledContours .* 255;
            window = ( (x>=loX & x<=hiX) & r<outerRad);
            tmpvar = zeros(m,n);
            tmpvar(window) = 1;
            tmpvar = repmat(tmpvar,[1 1 size(checks, 3)]);
            window = find(tmpvar == 1);
            
%             img = checks .* tmpvar;
%             img(~window) = bk;
            
            img         = bk*ones(size(checks));
            img(window) = checks(window);
        end
        % yet another loop to be able to move the checks...
        

        for ii=1:3
            images(:,:,(imgNum-1).*numMotSteps+ii)=uint8(img(:,:,ii));
        end
        for ii=numMotSteps-1:numMotSteps
            images(:,:,(imgNum-1).*numMotSteps+ii)=bk*ones(size(img(:,:,1)));
        end


        fprintf('.');drawnow;
    end
    fprintf('Done.\n');
end;



% Now we compile the actual sequence
% make stimulus sequence, make half and then add the rest as a flipped
% version of the first half
sequence = ones(duration.cycle.stimframes./(halfNumImages),1)*...
    [1:params.temporal.motionSteps:params.temporal.motionSteps*(halfNumImages)];

fixSeq=ones(duration.cycle.stimframes./(halfNumImages),1)*...
    [1:params.temporal.motionSteps:params.temporal.motionSteps*(halfNumImages)];
% fixSeqB=ones(duration.cycle.stimframes./2./halfNumImages,1);
% for ii=1:size(sequence, 2)
    motionSeq=ones(duration.cycle.stimframes./(halfNumImages.*2),1);
%     fixSeqB=ones(duration.cycle.stimframes./2./halfNumImages,1).*fixSeqB(end);

    motionSeq(numMotSteps-pauseDuration+1:numMotSteps)=numMotSteps;
    motionSeq(2*numMotSteps-pauseDuration+1:2*numMotSteps)=numMotSteps;
    motionSeq(3*numMotSteps-pauseDuration+1:3*numMotSteps)=numMotSteps;
    %motionSeq(4*numMotSteps-pauseDuration+1:4*numMotSteps)=numMotSteps;
%     if rand>oneBackFrequency
        motionSeq(numMotSteps+1:2*numMotSteps-pauseDuration)=2;
%         motionSeq(2*numMotSteps+1:3*numMotSteps-pauseDuration)=2;
%     else
%         if fixSeqB(numMotSteps)==1;
%             fixSeqB(numMotSteps+1:end)=2;
%         else
%             fixSeqB(numMotSteps+1:end)=1;
%         end
%     end
%     if rand>oneBackFrequency
        motionSeq(2*numMotSteps+1:3*numMotSteps-pauseDuration)=3;
        %motionSeq(3*numMotSteps+1:4*numMotSteps-pauseDuration)=4;
        %     else
%         if fixSeqB(2*numMotSteps)==1;
%             fixSeqB(2*numMotSteps+1:end)=2;
%         else
%             fixSeqB(2*numMotSteps+1:end)=1;
%         end
%     end
    motionSeq = repmat(motionSeq,1,size(sequence,2));
    %sequence(:)=sequence(:)+motionSeq-1;
    %fixSeq(:)=fixSeqB;
% end

sequence = sequence(:) + motionSeq(:) - 1;
fixSeq=fixSeq(:);

if params.insertBlanks.do,
    if params.insertBlanks.phaseLock, % keep 1 cycle and repeat
        completeCycle = sequence;
        sequence = [completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2))];
    else,
        sequence = repmat(sequence,params.ncycles,1);
        fixSeq = repmat(fixSeq,params.ncycles,1);
    end;
else,
    sequence = repmat(sequence,params.ncycles,1);
    fixSeq = repmat(fixSeq,params.ncycles,1);
end;

% we make only half so we need to flip the rest
% sep   = round(linspace(1,length(sequence)+1,5));
% rev = []; 
% fixrev=[];
% for n=1:4,
%     rev = [rev; flipud(sequence(sep(n):sep(n+1)-1))];
%     fixrev = [fixrev; flipud(fixSeq(sep(n):sep(n+1)-1))];
% end;
% sequence = [sequence; rev];
% fixSeq=[fixSeq; fixrev];


% direction        
if params.seqDirection~=0
	sequence = flipud(sequence);
end

% insert blanks (always of for 12 seconds)
% if params.insertBlanks.do,
%     seq2      = zeros(size(sequence));
%     oneCycle  = length(seq2)/params.insertBlanks.freq;
%     offTime   = ceil(10./params.tr).*params.tr; % make sure it's a multiple of the tr
%     offPeriod = ceil(offTime./duration.stimframe);
%     onPeriod  = oneCycle-offPeriod;
%     oneBlankCycle = [false(onPeriod,1); true(offPeriod,1)]; % blanks at the end
%     oneBlankCycle = [oneBlankCycle(oneCycle/2+1:end); oneBlankCycle(1:oneCycle/2)]; % blanks half way 
%     seq2      = repmat(oneBlankCycle,params.insertBlanks.freq,1);
%     add       = size(images,3)+1;
%     if isempty(params.loadMatrix),
%         for ii=2:length(seq2)
%             if seq2(ii)==1
%                 fixSeq(ii)=fixSeq(ii-1);
%             elseif seq2(ii)==0 && seq2(ii-1)==1 && fixSeq(ii)~=fixSeq(ii-1)
%                 fixSeq(ii:end)=3-fixSeq(ii:end);
%             end
%         end
%                 
%         
%         sequence(seq2)    = add;
%         images(:,:,add)   = uint8(ones(size(images,1),size(images,2)).*bk);
%     end;
%     clear seq2;
%     disp(sprintf('[%s]:Stimulus on for %.1f and off for %.1f seconds.',...
%         mfilename,onPeriod*duration.stimframe,offPeriod*duration.stimframe));
% end;
 if params.insertBlanks.do,
    seq2      = zeros(size(sequence));
    oneCycle  = length(seq2)/params.insertBlanks.freq;
    offTime   = ceil(10./params.tr).*params.tr; % make sure it's a multiple of the tr
    offPeriod = ceil(offTime./duration.stimframe);
    onPeriod  = oneCycle-offPeriod;
    oneBlankCycle = [false(onPeriod,1); true(offPeriod,1)]; % blanks at the end
    oneBlankCycle = [oneBlankCycle(oneCycle/2+1:end); oneBlankCycle(1:oneCycle/2)]; % blanks half way 
    seq2      = repmat(oneBlankCycle,params.insertBlanks.freq,1);
    add       = size(images,3)+1;
    if isempty(params.loadMatrix),
        sequence(seq2)    = add;
        images(:,:,add)   = uint8(ones(size(images,1),size(images,2)).*bk);
    end;
    clear seq2;
    disp(sprintf('[%s]:Stimulus on for %.1f and off for %.1f seconds.',...
        mfilename,onPeriod*duration.stimframe,offPeriod*duration.stimframe));
 end;


% fixation dot sequence
% change on the fastest every 6 seconds
minsec = 1.8./duration.stimframe;
fixSeq = ones(minsec,1)*round(rand(1,ceil(length(sequence)/minsec)));
fixSeq = fixSeq(:)+1;
fixSeq = fixSeq(1:length(sequence));
% % force binary
fixSeq(fixSeq>2)=2; 
fixSeq(fixSeq<1)=1;


% Insert the preappend images by copying some images from the
% end of the seq and tacking them on at the beginning
sequence = [sequence(length(sequence)+1-duration.prescan.stimframes:end); sequence];
timing   = [0:length(sequence)-1]'.*duration.stimframe;
cmap     = params.display.gammaTable;
fixSeq   = [fixSeq(length(fixSeq)+1-duration.prescan.stimframes:end); fixSeq];

% make stimulus structure for output
stimulus = createStimulusStruct(images,cmap,sequence,[],timing,fixSeq);

% save matrix if requested
if ~isempty(params.saveMatrix),
    save(params.saveMatrix,'images');
end;

