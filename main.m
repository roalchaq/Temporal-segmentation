
function main
   mainVideoTest
end


function mainVideoTest
    videoReader = VideoReader('Dissolve1-15 480x320 25 fps.mp4');   
    nFrames = videoReader.NumberOfFrames;
    frameRate = videoReader.FrameRate;
    videoWriter = VideoWriter('ADBF_Dissolve1-15 480x320 25 fps.avi', 'Uncompressed AVI');
    videoWriter.FrameRate = frameRate;
    open(videoWriter);
    for k = 1 : nFrames
            disp(['processing frame: ' num2str(k) ' of ' num2str(nFrames)]);            
            U = read(videoReader, k); 
            %to 360p
            %U = imresize(U, 0.5);
            %F1 = U;
            tic;
            F1 = processImage(U);   
            time = toc
            F1 = F1(20:size(F1,1)-20, 20:size(F1,2)-20, :);
            writeVideo(videoWriter, F1);
           %end
    end
    close(videoWriter);
    
end
%used parameters for the paper IWOBI 2014
function F = processImage(U)
    %asigne los parametros pertinentes, antes wRsize 21, sigma_r = 7
    wRSize = 21; sigma_s = wRSize/1.5; sigma_r = 10; lambda = 2;
    %fBilD = filterDeceivedBilateral(U, wRSize, sigma_s, sigma_r, lambda); 
    %fBilD = filterDGF(lambda, wRSize, 1, sigma_r, U);
    fMovAv = m2_filter(U,3,3);
    fBilD = uint8(fMovAv);

    F = fBilD;
end

%Input image must be from 0 to 255
function F = filterDeceivedBilateral(U, wSize, sigma_s, sigma_r, lambda)   
    %the image has to to have values from 0 to 1
    amps = [lambda*0.3 lambda lambda*0.2 ];  trap1 = [5 20 35 90];  trap2 = [70 100 150 255];  
   % amps = [lambda*0 lambda lambda ]; 
    
    %trap1 = [1 1 255 255];  trap2 = [254 254 255 255];  
    Unorm = double(U)/255;     
    [L, alfaMat, Vnorm] = adaptiveLaplacian(Unorm, amps, trap1, trap2);
    sigma = [sigma_s, sigma_r];
    F = bfilterDeceived(Unorm, L, wSize, sigma);   
    %putting back everything
    F = uint8(255 * F);
end

function [ fima] = filterDGF(lambda, wRSize, numIter, std_to_filter, U)
        
        amps = [lambda lambda lambda];  trap1 = [2 10 30 50];  trap2 = [30 100 150 255];  
		
        
        %ds = 0.01*diff(getrangefromclass(U)).^2;       
        ds = std_to_filter * std_to_filter;
        % Size of the image    
        tic;
            %in the cie lab color model
            L = adaptiveLaplacian(U, amps, trap1, trap2);
            G = colorspace('Lab<-RGB', U);
            %GF execution      
            fima = L;
            for i = 1 : numIter
                %preserves the original image as guidance
                %iterates to diminish noise
                %fimaAntes = uint8(fima);
                if(i >= 2)
                    wRSize = 2;
                end
                fima = imguidedfilter(fima, G, 'NeighborhoodSize',[wRSize*2 wRSize*2], 'DegreeOfSmoothing', ds);      
                fima = colorspace('RGB<-Lab', fima);
            end
        
        fima = uint8(255*fima);
end

%Shows the image on a window
function showOnWindow(U, name)
    figure; 
    T = uint8(U);
    imshow(T);
    title(name);
end


