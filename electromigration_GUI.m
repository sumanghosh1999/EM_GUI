function electromigration_GUI
    % Main GUI window for image selection
    f = figure('Name', 'Select an Interconnect Network Type', 'Position', [200 200 1200 500]);
    
    % Load images
    img1Data = imread('image1.png');
    img2Data = imread('image2.png');
    img1Data = imresize(imread('image1.png'), [300 500]);
    img2Data = imresize(imread('image2.png'), [300 500]);
    
    img1 = uicontrol('Style', 'pushbutton', 'CData', img1Data, 'Position', [50 100 500 300], 'Callback', @(~,~) setappdata(f, 'choice', 1));
    img2 = uicontrol('Style', 'pushbutton', 'CData', img2Data, 'Position', [600 100 500 300], 'Callback', @(~,~) setappdata(f, 'choice', 2));
    
    uicontrol('Style', 'pushbutton', 'String', 'Submit', 'Position', [525 20 100 30], 'Callback', @submitSelection);

    % Data storage
    global recordedData;
    recordedData = struct();

    % Wait for user interaction
    uiwait(f);

    function submitSelection(~,~)
        choice = getappdata(f, 'choice');
        if isempty(choice)
            errordlg('Please select an Interconnect Network Type');
        else
            recordedData.f = choice;
            close(f);
            if choice == 1
                inputValues1();
            else
                inputValues2();
            end
        end
    end

    function inputValues1
        prompt = {'m:', 'n:', 'I_injection (mA) [Magnitude]:', 'Temperature (K):'};
        title = 'Input parameters';
        dims = [3 60];
        definput = {'1', '1', '0', '300'};
        answer1 = inputdlg(prompt, title, dims, definput);
        if isempty(answer1), return; end
        m = str2double(answer1{1});
        n = str2double(answer1{2});
        I = str2double(answer1{3});
        T = str2double(answer1{4});

        metal = questdlg('Select Metal Type:', 'Metal Type', 'Cu', 'Al','Cu');
        if isempty(metal), return; end
        confinement = questdlg('Select Confinement Material:', 'Confinement Material', 'SiO2', 'Si3N4', 'SiO2');
        if isempty(confinement), return; end
        
        recordedData.m = m;
        recordedData.n = n;
        recordedData.I = I;
        recordedData.T = T;
        recordedData.metal = metal;
        recordedData.confinement = confinement;

        inputValues123(m, n);
    end

    function inputValues2
        prompt = {'n:', 'I_injection (mA) [Magnitude]:', 'Temperature (K):'};
        title = 'Input paramters';
        dims = [3 60];
        definput = {'1', '0', '300'};
        answer2 = inputdlg(prompt, title, dims, definput);
        if isempty(answer2), return; end
        n = str2double(answer2{1});
        I = str2double(answer2{2});
        T = str2double(answer2{3});

        metal = questdlg('Select Metal Type:', 'Metal Type', 'Cu', 'Al','Cu');
        if isempty(metal), return; end
        confinement = questdlg('Select Confinement Material:', 'Confinement Material', 'SiO2', 'Si3N4', 'SiO2');
        if isempty(confinement), return; end

        recordedData.n = n;
        recordedData.I = I;
        recordedData.T = T;
        recordedData.metal = metal;
        recordedData.confinement = confinement;

        inputValuesN(n);
    end

    function inputValues123(m, n)
        recordedData.segmentValues = inputMultipleValues(m, 'Left-Segment Parameters (L, W, H, j)');
        recordedData.centralValues = inputSingleValue('Central Segment Parameters (L, W, H, j)');
        recordedData.additionalSegmentValues = inputMultipleValues(n, 'Right-Segment Parameters (L, W, H, j)');
    end

    function inputValuesN(n)
        recordedData.segmentValues = inputMultipleValues(2^n-1, 'All Segment Parameters (L, W, H, j)');
    end

    function data = inputMultipleValues(numEntries, title)
        data = cell(numEntries, 4);
        for i = 1:numEntries
            prompt = {sprintf('Segment %d| L (um):', i), sprintf('Segment %d| W (um):', i), sprintf('Segment %d| H (um):', i), sprintf('Segment %d| j (GA/m^2):', i)};
            dims = [3 60];
            definput = {'100', '0.02366', '0.068', '1'};
            answer = inputdlg(prompt, title, dims, definput);
            if isempty(answer), return; end
            data(i, :) = answer;
        end
    end

    function data = inputSingleValue(title)
        prompt = {'L (um):', 'W (um):', 'H (um):', 'j (GA/m^2):'};
        dims = [3 60];
        definput = {'100', '0.02366', '0.068', '1'};
        data = inputdlg(prompt, title, dims, definput);
    end

end
