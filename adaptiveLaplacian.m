%Input: Image with values from 0 to 255
%Output: double image with undefined range and negatives
function [LAdaptive, alfaMat, Vnorm] = adaptiveLaplacian(U, amps, trap1, trap2)
    
    if(size(U, 3) == 1)
        U = double(U) ./ double(max(U(:)));
    else
        U = colorspace('Lab<-RGB', U);
       % U = double(U) ./ double(max(U(:)));
    end
    %Adaptive implementation of the laplacian filter
    %local variance
    V = stdfilt(U).^2;
    %Normalization without the borders
    S = 35;
    Vtemp = V(S : size( V, 1) - S, S : size(V, 2) - S);
    %OJO, ESTA NORMALIZADA RESPECTO AL MAXIMO DE V!
    Vnorm = uint8(255*(V./max(Vtemp(:))));
    [LAdaptive, alfaMat]= getAdaptiveLaplacian( Vnorm, U, amps, trap1, trap2);
   
end

%Calculates Adaptive laplacian image from a normalized local variance image
function [LAdaptive, alfaMat]= getAdaptiveLaplacian( Vnorm, U, amps, trap1, trap2) 
    F = generateTrapMembership(amps, trap1, trap2);
    Z = filterLaplacian2(U);   
    lim = max(U(:));
   % normalization
    maxZ = max(abs(Z(:)));
    Z = lim * (Z ./ maxZ);
    %unsharp masking
    %Laplacian of gaussian            
    LAdaptive = zeros(size(U(:,:,:))); 
    alfaMat = zeros(size(U(:,:, :))); 
    dimsU = size(U);
    if(length(dimsU) == 2)
        temp = ones(1,3);
        temp(1:2) = dimsU;
        dimsU = temp;
    end
    for i = 1: dimsU(1)
        for j = 1:dimsU(2)
            for k = 1:dimsU(3)
                varCurr = Vnorm(i, j, k);
                alfaCurr = getWeight(varCurr, F);
                alfaMat(i, j, k) = alfaCurr;
                %adaptive amplification
                LAdaptive(i, j, k) = U(i, j, k) + Z(i, j, k) * alfaCurr;
            end
        end
    end  
end
%applies a 11x11 laplacian mask (less sensitive to noise)
function F = filterLaplacian2(U)
    h =  fspecial('log', 17, 0.005);
    h =  -h;
    F = imfilter(U, h); 
end

%Generates the lambda laplacian gain functional from variance
function F = generateTrapMembership(amps, trap1, trap2)
    amp0 = amps(1);
    amp1 = amps(2);
    amp2 = amps(3);
    %from 0 to 255
    fTrap1 = trapmf(0:255, trap1);
    fTrap2 = trapmf(0:255, trap2);
    F = amp1 * fTrap1 + amp2 * fTrap2 + amp0;  
end



%Evaluates the variance to get the lambda coefficient
function alfa = getWeight(variance, F)
    variance = round(variance);
    if(variance >= 255)
        variance = 254;
    end

    %GMM = generateGMM(means, stds, maxValues, 255);    
    val =  variance + 1;
    alfa = F(val);
end




