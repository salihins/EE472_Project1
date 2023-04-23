function [Y_bus] = e230405_Aydin(cdf_path)
    dataFile = cdf_path;
    format long;
    [sBase] = readTitle(dataFile);
    [busStart,busEnd,branchStart,branchEnd] = getLineNum(dataFile);
    busNum = busEnd-busStart+1;
    busData = writeBusData(dataFile,busStart,busEnd);
    branchData = writeBranchData(dataFile,branchStart,branchEnd);
    Y_bus = writeOnBus(busData,branchData,busNum);
end

function [sBase] = readTitle(dataFile)
    rawData = readlines(dataFile);
    splitRawData = split(rawData(1));
    sBase = splitRawData(5);
end

function [numOfBus] = findYBusEntries(data,busNum,busData)  
    for i = 1:busNum
        if data == busData(i,1)
            numOfBus = i;
            return
        end
    end
end 
function [busStart,busEnd,branchStart,branchEnd] = getLineNum(dataFile)
    i = 1;
    busStart = 0;
    branchStart = 0;
    busEnd = 0;
    branchEnd = 0;
    rawData = readlines(dataFile);
    while i
        x = split(rawData(i));
        if x(1) == "BUS" && x(2) == "DATA"
            busStart = i+1;
        elseif x(1) == "-999" && branchStart == 0
            busEnd = i-1;            
        elseif x(1) == "BRANCH" && x(2) == "DATA"
            branchStart = i+1;          
        elseif x(1) == "-999" && branchStart ~= 0
            branchEnd = i-1;
            break
        end
        i = i+1;    
    end
end

function [busData] = writeBusData(dataFile,busStart,busEnd)
    rawData = readlines(dataFile);
    busData = {};
    for i = busStart:busEnd
        rowData = split(rawData(i));
        if busEnd-busStart+1 == 300
            if rowData(1) == ""
                busData = [busData; rowData(2),rowData(6),rowData(7),rowData(9),rowData(10),rowData(11),rowData(12),rowData(17),rowData(18)];
                busData = double(busData);
            else
                busData = [busData; rowData(1),rowData(5),rowData(6),rowData(8),rowData(9),rowData(10),rowData(11),rowData(16),rowData(17)];
                busData = double(busData);
            end

        elseif rowData(1) == "" && busEnd-busStart+1 ~= 300
            busData = [busData; rowData(2),rowData(end-13),rowData(end-12),rowData(end-10),rowData(end-9),rowData(end-8),rowData(end-7),rowData(end-2),rowData(end-1)];
            busData = double(busData);
        
        elseif rowData(1)~=""
            busData = [busData; rowData(1),rowData(end-13),rowData(end-12),rowData(end-10),rowData(end-9),rowData(end-8),rowData(end-7),rowData(end-2),rowData(end-1)];
            busData = double(busData);
        end
    end
end

function [branchData] = writeBranchData(dataFile,branchStart,branchEnd)
    rawData = readlines(dataFile);
    branchData = {};
    for i = branchStart:branchEnd
        rowData = split(rawData(i));
        if rowData(1) == ""
            branchData = [branchData; rowData(2),rowData(3),rowData(8),rowData(9),rowData(10),rowData(16),rowData(17)];
            branchData = double(branchData);
        else
            branchData = [branchData; rowData(1),rowData(2),rowData(7),rowData(8),rowData(9),rowData(15),rowData(16)];
            branchData = double(branchData);
        end
        
    end
end


function [yBus] = writeOnBus(busData,branchData,busNum)
    yBus = zeros(busNum,busNum);
    % Find Yij entries of Ybus matrix.
    for i=1:length(branchData)
        zBus = 0;
        zBus = branchData(i,3) + j * branchData(i,4);
        busEntity = -(1/zBus);
        if branchData(i,6) ~= 0
            busEntity = busEntity/branchData(i,6);
        end
        yBusX = findYBusEntries(branchData(i,1),busNum,busData);
        yBusY = findYBusEntries(branchData(i,2),busNum,busData);
        yBus(yBusX,yBusY) = yBus(yBusX,yBusY) + busEntity;
        yBus(yBusY,yBusX) = yBus(yBusY,yBusX) + busEntity;
    end

    % Fill the empty entries of lower triangle matrix.
%      yBus = yBus+transpose(yBus);



    % Sum all the columns to find Yii entry.
    for i=1:busNum
        temp = 0;
        for k= 1: busNum
            if i ~= k
                temp = temp + yBus(i,k);
            end
        end
        yBus(i,i) = yBus(i,i) - temp;
    end
    % Add B line to the Yii entries.
    for i=1:length(branchData)
        if branchData(i,5) ~= 0 && branchData(i,6) == 0
            zBus = 0;
            zBus = (j * branchData(i,5))/2;
            yBusX = findYBusEntries(branchData(i,1),busNum,busData);
            yBusY = findYBusEntries(branchData(i,2),busNum,busData);
            yBus(yBusX,yBusX) = yBus(yBusX,yBusX)+zBus;
            yBus(yBusY,yBusY) = yBus(yBusY,yBusY)+zBus;
        end
    end

    % Consider the tap ratios.
    for i=1:length(branchData)
        if branchData(i,6) ~= 0
            tapRatio = branchData(i,6);
            endCoeff = (tapRatio-1)/(tapRatio);
            fromCoeff = (1-tapRatio)/(tapRatio*tapRatio);
            zBus = 0;
            zBus = branchData(i,3) + j * branchData(i,4);
            busEntity = (1/zBus);
            addToFrom = busEntity.*fromCoeff;
            addToEnd = busEntity.*endCoeff;

            yBusX = findYBusEntries(branchData(i,1),busNum,busData);
            yBusY = findYBusEntries(branchData(i,2),busNum,busData);
            
            yBus(yBusX,yBusX) = yBus(yBusX,yBusX)+addToFrom;
            yBus(yBusY,yBusY) = yBus(yBusY,yBusY)+addToEnd;
        end
    end  


    % Add B shunt to the Yii entries
    for i=1:busNum
        if busData(i,8) ~= 0 || busData(i,9) ~= 0
            zBus = 0;
            zBus = busData(i,8) + j * busData(i,9);

            yBusX = findYBusEntries(busData(i,1),busNum,busData);
%             yBusY = findYBusEntries(busData(i,1),busNum,busData);
            yBus(yBusX,yBusX) = yBus(yBusX,yBusX)+zBus;
            
        end
    end
     
end
