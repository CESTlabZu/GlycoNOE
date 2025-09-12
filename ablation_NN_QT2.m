%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ABLATION II QUANTIFICATION NETWORK ARCHITECTURE
%
% Authors: Leqi Yin, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function net = ablation_NN_QT2(minLength, filterSize, numFilters)


    % Initialize empty layerGraph
    net = dlnetwork;

    inLayer = sequenceInputLayer(2, "MinLength", minLength, "Name", "input");
    net = addLayers(net, inLayer);






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



    ZS_dec = [

    % Decoder with skip connections
    transposedConv1dLayer(filterSize, numFilters, "Stride", 2, "Cropping", [1 1], "Name", "ZS_dec1")
    convolution1dLayer(filterSize, numFilters,"Padding", "same",  "Name", "ZS_decConv1")
    eluLayer("Name", "ZS_decElu1")

    transposedConv1dLayer(filterSize, numFilters/2, "Stride", 2, "Cropping", [1 1], "Name", "ZS_dec2")
    convolution1dLayer(filterSize, numFilters/2,"Padding", "same",  "Name", "ZS_decConv2")
    eluLayer("Name", "ZS_decElu2")

    transposedConv1dLayer(filterSize, numFilters/4, "Stride", 2, "Cropping", [1 1], "Name", "ZS_dec3")
    convolution1dLayer(filterSize, numFilters/4,"Padding", "same",  "Name", "ZS_decConv3")
    eluLayer("Name", "ZS_decElu3")

    convolution1dLayer(filterSize, 1, ...
    "Padding","same","Name","finalConv1")
    ];
    net = addLayers(net, ZS_dec);



    % LD Path (Lorentzian fitting) with consistent downsampling
    LD_dec = [

    transposedConv1dLayer(filterSize, numFilters, "Stride", 2, "Cropping", [1 1], "Name", "LD_dec1")
    convolution1dLayer(filterSize, numFilters,"Padding", "same",  "Name", "LD_decConv1")
    eluLayer("Name", "LD_decElu1")

    transposedConv1dLayer(filterSize, numFilters/2, "Stride", 2, "Cropping", [1 1], "Name", "LD_dec2")
    convolution1dLayer(filterSize, numFilters/2,"Padding", "same",  "Name", "LD_decConv2")
    eluLayer("Name", "LD_decElu2")

    transposedConv1dLayer(filterSize, numFilters/4, "Stride", 2, "Cropping", [1 1], "Name", "LD_dec3")
    convolution1dLayer(filterSize, numFilters/4,"Padding", "same",  "Name", "LD_decConv3")
    eluLayer("Name", "LD_decElu3")

    convolution1dLayer(filterSize, 1, ...
    "Padding","same","Name","finalConv2")


    ];
    net = addLayers(net, LD_dec);


    final_conv = [
        depthConcatenationLayer(2, "Name", "finalConcat");
        convolution1dLayer(filterSize, 1, ...
        "Padding","same","Name","finalConv3")
        ];

    net = addLayers(net, final_conv);


    % Connect Layers
    net = connectLayers(net, "input", "shared_enc1");


    net = connectLayers(net, "bottleneck", "ZS_dec1");
    net = connectLayers(net, "bottleneck", "LD_dec1");
    net = connectLayers(net, "finalConv1","finalConcat/in1");
    net = connectLayers(net, "finalConv2","finalConcat/in2");

end
