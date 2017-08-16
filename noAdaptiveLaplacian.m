function L = noAdaptiveLaplacian(U, lambda)
    U = colorspace('Lab<-RGB', U); 
    L = filterUM_laplacianLAB(U, lambda);
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