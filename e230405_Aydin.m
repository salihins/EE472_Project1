function [V_bus,Angle_bus,Pg, Qg] = e230405_Aydin(cdf_path)
    dataFile = cdf_path;
    format long;
    epsilon = 0.001;
    [busStart,busEnd,branchStart,branchEnd] = getLineNum(dataFile);
    busData = writeBusData(dataFile,busStart,busEnd);    
    sBase = readTitle(dataFile);
    sBase = str2num(sBase);
    busNum = busEnd-busStart+1;
    branchData = writeBranchData(dataFile,branchStart,branchEnd);
    Y_bus = writeOnBus(busData,branchData,busNum);
    slackBus = findSlack(busData);
    
    slackBusNumber = findYBusEntries(slackBus,busNum,busData);
    pvBuses = findPVBus(busData);
    pqBuses = findPQBus(busData);
    
    Y_bus = writeOnBus(busData,branchData,busNum);
    bMatrix = -imag(Y_bus);
    bPrime = createBprime(Y_bus,slackBusNumber);
    bDoublePrime = createBDoublePrime(Y_bus,pqBuses,pvBuses,slackBus,bPrime,busNum,busData);
    busVoltage = findBusVoltages(busData,slackBusNumber,busNum);
    for t=1:length(pqBuses)
        busVoltage(pqBuses(t)) = 1;
    end
    [pSpec,qSpec] = findSpecValues(busData,slackBusNumber,sBase);
    initA = zeros(busNum,1);
%     tic;
    for i = 1:20
        [deltaP,P,deltaQ,Q] = findDeltaPQ(Y_bus,pSpec,qSpec,busVoltage,initA);
        firstCorMat = deltaP./busVoltage;
        firstCorMat(slackBusNumber,:) = [];
        deltaA = inv(bPrime)*(firstCorMat);
        secondCorMat = deltaQ./busVoltage;
        thirdCorMat = zeros(length(pqBuses),1);
        for k=1:length(pqBuses)
            thirdCorMat(k) = secondCorMat(pqBuses(k),:);
        end
        deltaV = inv(bDoublePrime)*(thirdCorMat);
        

        if abs(deltaP(2:end))<=epsilon
            fprintf("deltaP Converges\n")
            if sum(abs(deltaQ(2:end))<=epsilon)
                fprintf("deltaQ also converges\n")
                V_bus = busVoltage;
                Angle_bus = initA * (360/(2*pi));
                Pg = P;
                Qg = Q;

                for t=1:length(pqBuses)
                    Pg(pqBuses(t)) = 0;
                    Qg(pqBuses(t)) = 0;
                end
             break
             end
        end
        fprintf("iteration:%d\n",i)
        for p=1:(busNum-1)
            initA(p+1) = initA(p+1)+deltaA(p);
        end
               
        for t=1:length(pqBuses)
            busVoltage(pqBuses(t)) = busVoltage(pqBuses(t)) + deltaV(t);
        end

    end
%     toc;
    Pgen = Pg(1);
    Qgen = Qg(1);
    Pg = zeros(length(deltaP),1);
    Qg = zeros(length(deltaQ),1);
    Pg(1) = P(1);
    Qg(1) = Q(1);
    for i = 1:length(pvBuses)
        if busData(pvBuses(i),4)>0
            Pgen = Pgen+ P(pvBuses(i)) + busData(pvBuses(i),4)/sBase;
            Pg(pvBuses(i)) = P(pvBuses(i)) +busData(pvBuses(i),4)/sBase;
        else
            Pgen = Pgen+ P(pvBuses(i));
            Pg(pvBuses(i)) = P(pvBuses(i));
        end

        if busData(pvBuses(i),5)>0
          Qgen = Qgen+ Q(pvBuses(i)) +busData(pvBuses(i),5)/sBase;
          Qg(pvBuses(i)) = Q(pvBuses(i)) +busData(pvBuses(i),5)/sBase;
        else
          Qgen = Qgen+ Q(pvBuses(i));
          Qg(pvBuses(i)) = Q(pvBuses(i));
        end
    end
  Pload = sum(busData(1:end,4))/sBase;
  Qload = sum(busData(1:end,5))/sBase;
  PowerLoss = Pgen-Pload;
  fprintf("Total Loss is %d MW",PowerLoss*sBase)
  

end

function [deltaP,P,deltaQ,Q] = findDeltaPQ(Y_bus, initP,initQ,initV,initA)
    P = zeros(length(Y_bus),1);
    Q = zeros(length(Y_bus),1);
    G = real(Y_bus);
    B = imag(Y_bus);
    for i = 1:length(Y_bus)
        for k = 1:length(Y_bus)
            P(i) = P(i) + initV(k)*initV(i)*((G(i,k)*cos(initA(i)-initA(k)))+(B(i,k)*sin(initA(i)-initA(k))));
        end
    end
    deltaP = initP-P;
        
    for i = 1:length(Y_bus)
        for k = 1:length(Y_bus)
            Q(i) = Q(i) + initV(k)*initV(i)*((G(i,k)*sin(initA(i)-initA(k)))-(B(i,k)*cos(initA(i)-initA(k))));
        end
    end

    deltaQ = initQ-Q;
    
end

function [pSpec,qSpec]= findSpecValues(busData,slackBusNumber,sBase)

    pSpec = (busData(1:end,6)-busData(1:end,4))/sBase;
    qSpec = (busData(1:end,7)-busData(1:end,5))/sBase;
%     pSpec(slackBusNumber,;) = [];
%     qSpec(slac
end
function [V] = findBusVoltages(busData,slackBusNumber,busNum)
    
    V = ones(busNum,1);
    for i=1:length(busData)
        V(i,1) = busData(i,3);
    end
    
%      V(slackBusNumber,:) = [];
    
end
function pvBus = findPVBus(busData)
    pvBus = [];
    for i = 1:length(busData)
        if busData(i,2) == 2
            pvBus = [pvBus,busData(i,1)];
        end
    end
end

function pqBus = findPQBus(busData)
    pqBus = [];
    for i = 1:length(busData)
        if busData(i,2) == 0
            pqBus = [pqBus,busData(i,1)];
        end
    end
end

function [sBase] = readTitle(dataFile)
    rawData = readlines(dataFile);
    splitRawData = split(rawData(1));
    sBase = splitRawData(5);
end

function slackBusNum = findSlack(busData)
    for i = 1:length(busData)
        if busData(i,2) == 3
            slackBusNum = busData(i,1);
            break;
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
function [numOfBus] = findYBusEntries(data,busNum,busData)  
    for i = 1:busNum
        if data == busData(i,1)
            numOfBus = i;
            return
        end
    end
end

function Bprime = createBprime(Y_bus,slackBusNumber)
    B = -imag(Y_bus);
    Bcopy = B;
    
    Bcopy(slackBusNumber,:) = [];
    Bcopy(:,slackBusNumber) = [];
    
    Bprime = Bcopy;

end


function BDoubleprime = createBDoublePrime(Y_bus, pqBuses,pvBuses,slackBus,bPrime,busNum,busData)
      bPrimeCopy = bPrime;
      pqBusNumber = length(pqBuses);
      pvBusNumber = length(pvBuses);
      slackBusNumber = findYBusEntries(slackBus,busNum,busData);
      deleteCount = 0;
      for i=1:pvBusNumber
          numOfBus = findYBusEntries(pvBuses(i),busNum,busData);
          
          if numOfBus < slackBusNumber
              bPrimeCopy(numOfBus-deleteCount,:) = [];
              bPrimeCopy(:,numOfBus-deleteCount) = [];
          else
              bPrimeCopy(numOfBus-1-deleteCount,:) = [];
              bPrimeCopy(:,numOfBus-1-deleteCount) = [];
          end
          deleteCount = deleteCount + 1;
      end
      
      BDoubleprime = bPrimeCopy;
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
