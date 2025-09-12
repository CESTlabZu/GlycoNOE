%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GLARES QUANTIFICATION NETWORK ARCHITECTURE
%
% Authors: Leqi Yin, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function net = GLARES_NN_QT(minLength, filterSize, numFilters)


    % Initialize empty layerGraph
    net = dlnetwork;

    % Input layer: 2 channels
    inLayer = sequenceInputLayer(2, "MinLength", minLength, "Name", "input");
    net = addLayers(net, inLayer);

    % split the 2 channels
    splitLayer1 = SplitChannelsLayer3("split1");
    net = addLayers(net, splitLayer1);









    shared_enc = [
        % Encoder
        convolution1dLayer(filterSize, numFilters/4, "Padding", "same", "Name", "shared_enc1")
        eluLayer("Name", "shared_encElu1")
        maxPooling1dLayer(2,'Stride',2,"Padding","same")

        convolution1dLayer(filterSize, numFilters/2, "Padding", "same",  "Name", "shared_enc2")
        eluLayer("Name", "shared_encElu2")
        maxPooling1dLayer(2,'Stride',2,"Padding","same")

        convolution1dLayer(filterSize, numFilters, "Padding", "same",  "Name", "shared_enc3")
        eluLayer("Name", "shared_encElu3")
        maxPooling1dLayer(2,'Stride',2,"Padding","same")

        flattenLayer("Name", "flatten")
        fullyConnectedLayer(numFilters/2, "Name", "bottleneck")
        ];

    net = addLayers(net, shared_enc);



    main_dec = [

    % Decoder with skip connections
    transposedConv1dLayer(filterSize, numFilters/4, "Stride", 2, "Cropping", [1 1], "Name", "main_dec1")
    convolution1dLayer(filterSize, numFilters/4,"Padding", "same",  "Name", "main_decConv1")
    eluLayer("Name", "main_decElu1")

    transposedConv1dLayer(filterSize, numFilters/8, "Stride", 2, "Cropping", [1 1], "Name", "main_dec2")
    convolution1dLayer(filterSize, numFilters/8,"Padding", "same",  "Name", "main_decConv2")
    eluLayer("Name", "main_decElu2")

    transposedConv1dLayer(filterSize, numFilters/16, "Stride", 2, "Cropping", [1 1], "Name", "main_dec3")
    convolution1dLayer(filterSize, numFilters/16,"Padding", "same",  "Name", "main_decConv3")
    eluLayer("Name", "main_decElu3")

    convolution1dLayer(filterSize, 1, ...
    "Padding","same","Name","finalConv1")
    ];
    net = addLayers(net, main_dec);



    % Aux Path (Lorentzian fitting) with consistent downsampling
    Aux_dec = [

    % Decoder with skip connections
    transposedConv1dLayer(filterSize, numFilters/4, "Stride", 2, "Cropping", [1 1], "Name", "Aux_dec1")
    convolution1dLayer(filterSize, numFilters/4,"Padding", "same",  "Name", "Aux_decConv1")
    eluLayer("Name", "Aux_decElu1")

    transposedConv1dLayer(filterSize, numFilters/8, "Stride", 2, "Cropping", [1 1], "Name", "Aux_dec2")
    convolution1dLayer(filterSize, numFilters/8,"Padding", "same",  "Name", "Aux_decConv2")
    eluLayer("Name", "Aux_decElu2")

    transposedConv1dLayer(filterSize, numFilters/16, "Stride", 2, "Cropping", [1 1], "Name", "Aux_dec3")
    convolution1dLayer(filterSize, numFilters/16,"Padding", "same",  "Name", "Aux_decConv3")
    eluLayer("Name", "Aux_decElu3")

    convolution1dLayer(filterSize, 1, ...
    "Padding","same","Name","finalConv2")

    depthConcatenationLayer(2, "Name", "mergeConcat1");

    convolution1dLayer(2, 1, ...
    "Padding","same","Name","Aux_Conv")
    ];
    net = addLayers(net, Aux_dec);


    finalconcatLayer = depthConcatenationLayer(2, "Name", "finalConcat");
    net = addLayers(net, finalconcatLayer);


    % Connect Layers
    net = connectLayers(net, "input", "split1");

    % main Path Connections
    net = connectLayers(net, "split1/ch2", "mergeConcat1/in2");
    net = connectLayers(net, "split1/mainch", "shared_enc1");
    net = connectLayers(net, "bottleneck", "main_dec1");
    net = connectLayers(net, "bottleneck", "Aux_dec1");
    net = connectLayers(net, "finalConv1","finalConcat/in1");
    net = connectLayers(net, "Aux_Conv","finalConcat/in2");

end
