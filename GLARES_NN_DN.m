%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GLARES DENOISING NETWORK ARCHITECTURE
%
% Authors: Leqi Yin, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function net = GLARES_NN_DN(minLength, filterSize, numFilters)


    % Initialize empty layerGraph
    net = dlnetwork;

    inLayer = sequenceInputLayer(5, "MinLength",minLength, "Name", "input");
    net = addLayers(net, inLayer);


    autoenc = [
        % Encoder
        convolution1dLayer(filterSize, numFilters/8, "Padding", "same", "Name", "conv1")
        eluLayer("Name", "encElu1")
        maxPooling1dLayer(2,'Stride',2,"Padding","same")

        convolution1dLayer(filterSize, numFilters/4, "Padding", "same",  "Name", "conv2")
        eluLayer("Name", "encElu2")
        maxPooling1dLayer(2,'Stride',2,"Padding","same")

        convolution1dLayer(filterSize, numFilters,"Padding", "same",  "Name", "conv3")
        eluLayer("Name", "encElu3")
        maxPooling1dLayer(2,'Stride',2,"Padding","same")

        flattenLayer("Name", "flatten")
        fullyConnectedLayer(numFilters/4, "Name", "bottleneck")

        transposedConv1dLayer(filterSize, numFilters, "Stride", 2, "Cropping", [1 1], "Name", "dec1")
        convolution1dLayer(filterSize, numFilters,"Padding", "same",  "Name", "conv4")
        eluLayer("Name", "decElu1")
    

        transposedConv1dLayer(filterSize, numFilters/4, "Stride", 2, "Cropping", [1 1], "Name", "dec2")
        convolution1dLayer(filterSize, numFilters/4,"Padding", "same",  "Name", "conv5")
        eluLayer("Name", "decElu2")
  

        transposedConv1dLayer(filterSize, numFilters/8, "Stride", 2, "Cropping", [1 1], "Name", "dec3")
        convolution1dLayer(filterSize, numFilters/8,"Padding", "same",  "Name", "conv6")
        eluLayer("Name", "decElu3")

        convolution1dLayer(filterSize, 5, "Padding","same","Name","finalConv1")



        ];

    net = addLayers(net, autoenc);




    % Connect Layers
    net = connectLayers(net, "input", "conv1");

