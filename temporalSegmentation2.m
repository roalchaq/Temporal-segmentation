

function temporalSegmentation2
   global distanceFunctions;
    j = 1;
    %Specify the distance functions to use
    distanceFunctions{j} = 'jeffrey_divergence'; j = j + 1;
    distanceFunctions{j} = 'jensen_shannon_divergence'; j = j + 1;
    distanceFunctions{j} = 'kolmogorov_smirnov_distance'; j = j + 1;
    distanceFunctions{j} = 'kullback_leibler_divergence'; j = j + 1;
    distanceFunctions{j} = 'quadratic_form_distance'; j = j + 1;
    distanceFunctions{j} = 'histogram_intersection'; j = j + 1;
    distanceFunctions{j} = 'chi_square_statistics'; j = j + 1;
    %execute the test
     folder = '';    
     [GT_Vector, gtTransitions]= readGroundtruth('anotaciones video sintetico.txt',25);
     psnrs(1, :) = compareBhattHistsVideos([folder 'Dissolve1-15_ADBF.avi'], [folder 'Dissolve1-15.mp4'], gtTransitions);
end




%groundtruth, a matrix with each row containing the cut frames
function [psnrs] = compareBhattHistsVideos(videoNameDBF, videoNameNoDBF, groundtruth)
    [distCoefsNoDBF, fps]  = getModBhattacharyyaDistFrames(videoNameNoDBF);
    [distCoefsDBF, fps]  = getModBhattacharyyaDistFrames(videoNameDBF);
    for i = 1 : size(distCoefsNoDBF, 2)
        bhattacharyyaCoefsDBF = distCoefsDBF(:, i);
        bhattacharyyaCoefsNoDBF = distCoefsNoDBF(:, i);
        
        thresholdedBhattAllDBF = thresholdBhattacharyya(bhattacharyyaCoefsDBF, 0, fps);
        [fpDBF(i), fnDBF(i)] = measureFpsAndFns(groundtruth, thresholdedBhattAllDBF);

        
        thresholdedBhattAllNoDBF = thresholdBhattacharyya(bhattacharyyaCoefsNoDBF, 0, fps);
        [fpNoDBF(i), fnNoDBF(i)] = measureFpsAndFns(groundtruth, thresholdedBhattAllNoDBF);
        
        
        psnrDBF = getPSNR(bhattacharyyaCoefsDBF, groundtruth);
        psnrNoDBF = getPSNR(bhattacharyyaCoefsNoDBF, groundtruth);
        
        psnrs(i, 1) = psnrNoDBF;
        psnrs(i, 2) = psnrDBF;
        temp = psnrDBF - psnrNoDBF;
        psnrs(i, 3) = temp / psnrNoDBF;
    end
    save('distanceResults');
    
end

function psnr = getPSNR(bhattacharyyaCoefs, groundtruth)
    maxSignal = 0;
    bhattacharyyaCoefs = bhattacharyyaCoefs * 100;
    bhattacharyyaCoefsNoise = bhattacharyyaCoefs;
    for i = 1:size(groundtruth, 1)
        currSignalInter = groundtruth(i, :);
        maxSignal = max([maxSignal  max(bhattacharyyaCoefs(currSignalInter(1):currSignalInter(2)))]);
        bhattacharyyaCoefsNoise(currSignalInter(1):currSignalInter(2)) = 0;        
    end
    %in no transition regions, signal should be zero
    meanNoise = mean(bhattacharyyaCoefsNoise(bhattacharyyaCoefsNoise~=0).^2);
    psnr = 20 * log10(maxSignal / sqrt(meanNoise));
end

%Favor implementar este para leer el archivo generado por el programa de
%andres
%en este caso es un archivo hard codeado
%se necesita que lea un archivo especifico y pase el momento del cut o
%dissolve en segundos
%debe retornar el tiempo en segundos del cut
%ademas del tipo de cut

function [GT_Vector, gtTransitions]= readGroundtruth(fileNameGT,fps)
fd = fopen(fileNameGT);
Data_raw = textscan(fd,'%*s%s%*s','Delimiter',':,');
fclose(fd);
deltaCut = 5;
deltaShort = 20;
deltaLong = 20;
Object = Class;

limite = length(Data_raw{1})-4;
m=1;
    for k = 3:6:limite
            Object.InitialFrame = str2double(Data_raw{1}{k});
            Object.LastFrame = str2double(Data_raw{1}{k+1});
            Object.ObjectType = Data_raw{1}{k+2};
            strType = Object.ObjectType;
            initialF = -1;
            if( strcmp(strType, '"Long Dissolve"'))
                initialF = Object.InitialFrame - deltaLong;
                last = Object.LastFrame + deltaLong;
            end
            if( strcmp(strType, '"Short Dissolve"'))
                initialF = Object.InitialFrame - deltaShort;
                last = Object.LastFrame + deltaShort;
            end
            if( strcmp(strType, '"Cut"'))
                initialF = Object.InitialFrame - deltaCut;
                last = Object.LastFrame + deltaCut;
            end
            %gtCutSeconds(m) = Object.InitialFrame*fps;
            if(initialF ~= -1)
                gtTransitions(m, :) = [initialF last];
                Object.SceneType = Data_raw{1}{k+3};
                GT_Vector(m) = Object;
                m=m+1;
            end

    end
end

%timeThresh, maximum tolerated 
%groundtruth cut positions in seconds, manually generated
function [fp, fn, tp] = measureFpsAndFns(gt, thresholdedBhattAll)
    % for each found cut, find 
    fn = 0;
    tp = 0;
   
    foundPositions = zeros(size(gt));
    for i = 1: length(gt)
        minGT = gt(i, 1);
        maxGT = gt(i, 2);
        if(sum(thresholdedBhattAll(minGT : maxGT))> 0)
          thresholdedBhattAll(minGT : maxGT) = 0;
          tp = tp + 1;
        else
           %thresholdedBhattAll(i)
           fn = fn + 1;
        end
        
    end    
    fp = sum(thresholdedBhattAll); 
    %plot(thresholdedBhattAll);
end


function [distCoefs, fps] = getModBhattacharyyaDistFrames(videoName)
    global distanceFunctions;
    videoReader = VideoReader(videoName);    
    nFrames = videoReader.NumberOfFrames;    
    fps = videoReader.frameRate;
    histH = 0;  
    histV = 0;
    j = 1;
    for k = 1 : nFrames
        disp(['processing frame: ' num2str(k) ' of ' num2str(nFrames)]);            
        U = read(videoReader, k); 
        %to 360p
        %U = imresize(U, 0.25);
        %imshow(U);
        nBins = 256; 
        normalize = 1;
        Uhsv = rgb2hsv(U);
        Uh = Uhsv(:,:,1);
        Uh = uint16(255 * Uh);
        Uv = uint16(255*Uhsv(:,:,3));
        Uvar = stdfilt(Uv);
        Uvar = uint8(255 * (Uvar ./ max(Uvar(:))));
        mask = ones(size(U));
        prevHistH = histH;
        prevHistV = histV;
        
        histH = Hhist(Uh, mask, nBins, normalize);
        histV = Hhist(Uvar, mask, nBins, normalize);
        
        if(sum(prevHistV) > 0)
            %implement here the custom bhattacharyya distance
            bhattacharyyaH = bhattacharyya(prevHistH, histH);
            bhattacharyyaV = bhattacharyya(prevHistV, histV);
            distCoefs(j, 1) = bhattacharyyaV * bhattacharyyaH;
            %iterates through every distance function to calculate it
            for i = 2:length(distanceFunctions)
                distCoefs(j, i) =  calculateDistance(prevHistH', histH', prevHistV', histV', distanceFunctions{i});
            end
            j = j + 1;
        end      
    end    
end

function dist = calculateDistance(prevHistH, histH, prevHistV, histV, funcDistName)   
    prevHistH = prevHistH + eps;
    histH = histH + eps;
    prevHistV = prevHistV + eps;
    histV = histV + eps;
    funcDist = str2func(funcDistName);
    if(strcmp(funcDistName, 'quadratic_form_distance'))
        distH = funcDist(prevHistH, histH, 1);
        distV = funcDist(prevHistV, histV, 1);
    else
        distH = funcDist(prevHistH, histH);
        distV = funcDist(prevHistV, histV);
    end
    dist = distH * distV;
end




function dist = emdHists(prevHistH, histH, prevHistV, histV)
    [signPrevH, valsPrevH] = histogram2Signature2(prevHistH);
    [signH, valsH] = histogram2Signature2(histH);
    [signPrevV, valsPrevV] = histogram2Signature2(prevHistV);
    [signV, valsV] = histogram2Signature2(histV);
    [x fvalH] = emd(signPrevH, signH, valsPrevH, valsH, @gdf);
    [x fvalV] = emd(signPrevV, signV, valsPrevV, valsV, @gdf);
    dist = fvalH * fvalV;
end



function [sign, vals] = histogram2Signature2(hist)
    j = 1;
    for i = 1:length(hist)
        if(hist(i)>0)
          sign(j) = i;
          vals(j) = hist(j);
          j = j + 1;
        end
    end 
    sign = sign';
    vals = vals';
end

function thresholdedBhattAll = thresholdBhattacharyya(bhattacharyyaCoefs, deltaMinutLocalThresh, fps)
    deltaFrames = length(bhattacharyyaCoefs) - 1;
    if(deltaMinutLocalThresh > 0)
        deltaFrames = deltaMinutLocalThresh * 60 * fps;
    end
    firstFrame = 1;
    finalFrame = firstFrame + deltaFrames;
    numDeltas = round(length(bhattacharyyaCoefs)/deltaFrames);
    thresholdedBhattAll = zeros(size(bhattacharyyaCoefs));
    for i = 1 : numDeltas
        bhattacharyyaCoefsLocal = bhattacharyyaCoefs(firstFrame:finalFrame);
        
        thresholdedBhatt = calculateThreshold(bhattacharyyaCoefsLocal);
        thresholdedBhattAll(firstFrame:finalFrame) = thresholdedBhatt;
        firstFrame = finalFrame;
        finalFrame = firstFrame + deltaFrames + 1;
    end
    if(finalFrame < length(bhattacharyyaCoefs) && deltaMinutLocalThresh > 0)
        bhattacharyyaCoefsLocal = bhattacharyyaCoefs(firstFrame:length(bhattacharyyaCoefs));
        thresholdedBhatt = calculateThreshold(bhattacharyyaCoefsLocal);
        thresholdedBhattAll(firstFrame:length(bhattacharyyaCoefs)) = thresholdedBhatt;
    end

end

function thresholdedBhatt = calculateThreshold(bhattacharyyaCoefsLocal)
    meanBhattCoefs = mean(bhattacharyyaCoefsLocal);
    stdBhattCoefs = std(bhattacharyyaCoefsLocal);
    td = meanBhattCoefs + stdBhattCoefs;
    tc = 2 * td;
    thresholdedBhatt = zeros(size(bhattacharyyaCoefsLocal));
    thresholdedBhatt(bhattacharyyaCoefsLocal > tc) = 2;
    thresholdedBhatt(bhattacharyyaCoefsLocal < tc & bhattacharyyaCoefsLocal > td) = 1;
    
end

%0 to 255 values to build histogram in input image
function H = Hhist(I, mask, nBins,Nind)
      

    if(nargin<3)
        Nind=0;
        % Default is un-normalized histogram
    end
    H=zeros(nBins, 1);
    for i=1:size(I,1)
        for j=1:size(I,2)
            if(mask(i, j) == 1)   
                val = I(i,j) + 1;               
                H(val) = H(val) + 1; 
            end
        end
    end
    % Un-Normalized histogram
    if(Nind==1)
        H = H./sum(H);
        % l1 normalization    
    end
end


function bdist = bhattacharyya(histogram1H, histogram2H)
    bins = size(histogram1H, 1);
    mean1H = double(mean(histogram1H));
    mean2H = double(mean(histogram2H));
    bcoeffSum = 0;
    for i=1:bins
        bcoeffSum = bcoeffSum + sqrt(double(histogram1H(i) * histogram2H(i)));
    end
    c1 = 1.0 / (sqrt(mean1H * mean2H * bins * bins));
    bdist = sqrt(1.0 - c1 * bcoeffSum);
    if(bdist == Inf)
        %bcoeff = 1;
        disp('suave');
    end

end