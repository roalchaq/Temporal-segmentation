
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-process input and select appropriate filter.
function B = bfilterDeceived(A, L, w, sigma)  
    if ~isfloat(A) || ~sum([1,3] == size(A,3)) || ...
          min(A(:)) < 0 || max(A(:)) > 1
       error(['Input image A must be a double precision ',...
              'matrix of size NxMx1 or NxMx3 on the closed ',...
              'interval [0,1].']);      
    end
    % Apply either grayscale or color bilateral filtering.
    if size(A,3) == 3
       B = bfltColorDeceived(A, L, w, sigma(1),sigma(2));     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements bilateral filter for color images.
%sigma range is multiplied by 100
function B = bfltColorDeceived(A, L, w, sigma_d,sigma_r)    
    % Convert input sRGB image to CIELab color space. 
    disp('using the CIELAB color space')
    %CIE LAB AB values go from 0 to 100
    A = colorspace('Lab<-RGB', A); 
   % L =  filterUM_laplacianLAB(A, lambda);  
    B = bfil2LAB_deceived(A, L, w, sigma_d,sigma_r);
    % Convert filtered image back to sRGB color space.   
    B = colorspace('RGB<-Lab',B);
end



function B = bfil2LAB_deceived(A, Laplacian, w, sigma_d,sigma_r)
    % Pre-compute Gaussian domain weights.
    [X,Y] = meshgrid(-w:w,-w:w);
    G = exp(-(X.^2 + Y.^2)/(2 * sigma_d^2));

    % Rescale range variance (using maximum luminance).
    %sigma range is multiplied by 100
    
    % Create waitbar.
    h = waitbar(0,'Applying the deceived bilateral filter...');
    set(h,'Name','Deceived Bilateral Filter Progress');

    % Apply bilateral filter.
    dim = size(A);
    B = zeros(dim);
   
    for i = 1:dim(1)
       for j = 1:dim(2)
             % Extract local region.
             iMin = max(i - w,1);
             iMax = min(i + w,dim(1));
             jMin = max(j - w,1);
             jMax = min(j + w,dim(2));
             I = A(iMin:iMax,jMin:jMax,:);
             %Laplacian window
             L = Laplacian(iMin:iMax,jMin:jMax,:);     
             % Compute Gaussian range weights.
             %done in the three layers
             dL = I(:,:,1) - A(i,j,1);
             da = I(:,:,2) - A(i,j,2);
             db = I(:,:,3) - A(i,j,3);
             H = exp(-(dL.^2 + da.^2 + db.^2) / (2 * sigma_r^2));

             % Calculate bilateral filter response.
             F = H .* G((iMin:iMax) - i + w + 1,(jMin:jMax) -j + w + 1);
             norm_F = sum(F(:));
             %The laplacian deceive consists on weighting the gaussian
             %function with the original image, and using the image values
             %of the laplacian image.
             B(i,j,1) = sum(sum(F.* L(:,:,1)))/norm_F;
             B(i,j,2) = sum(sum(F.* L(:,:,2)))/norm_F;
             B(i,j,3) = sum(sum(F.* L(:,:,3)))/norm_F;

       end
       waitbar(i/dim(1));
    end
    % Close waitbar.
    close(h);
end

%Regular non adaptive unsharped masking with laplacian filter
%normalizacion incorrecta
%Tercer metodo a utilizar
function F = filterUM_laplacianLAB(U, lambda1)
    Z = filterLaplacian2(U);   
    lim = max(U(:));
   % normalization
    maxZ = max(abs(Z(:)));
    Z = lim * (Z ./ maxZ);
    %unsharp masking
    F = U + lambda1 .* Z; 

end
%applies a 11x11 laplacian mask
function F = filterLaplacian2(U)
    h =  fspecial('log', 17, 0.005);
    h =  -h;
    F = imfilter(U, h); 
end




