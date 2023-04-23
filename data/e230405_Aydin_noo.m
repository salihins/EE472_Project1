function [V,Angle,Pgen_matrix,Qgen_matrix] = e230405_Aydin_old(cdf_path)
    epsilon = 0.001;
    dataFile = cdf_path;
    format long;
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
%     bPrime = createBprime(Y_bus,slackBusNumber);
%     bDoublePrime = createBDoublePrime(Y_bus,pqBuses,pvBuses,slackBus,bPrime,busNum,busData);
    [B,bPrime,bDoublePrime] = findB(Y_bus,slackBusNumber,pvBuses,pqBuses);
    x = 1;
    [Pspec,Qspec]=  specValues(busData,slackBusNumber,sBase);


  AngleData = zeros(length(Y_bus),1);
  V = busData(1:end,3);
  for i = 1:length(pqBuses)
      V(pqBuses(i)) = 1;
  end
  Vold = [V];

  Angle_old =[];
  deltaP_old =[];
  deltaQ_old = [];
  k=1;
  while k<21

    [deltaP,P] = findDeltaP(Y_bus,slackBusNumber,Pspec,sBase,AngleData,V);

      if sum(abs(deltaP(2:end))<epsilon) == length(deltaP)-1
          fprintf("deltaP Converges\n")
          if sum(abs(deltaQ_compare)<epsilon) == length(deltaQ_compare)
              fprintf("deltaQ also converges\n")
              k = 999;
          end
      else
          AngleData = FDLF_Angle(bPrime,bDoublePrime,deltaP,AngleData,V,slackBusNumber,pqBuses);
          Angle_old = [Angle_old AngleData];
      end
      %%%%%%%%%%%%%%%%%%%%%%
      [deltaQ,Q] = findDeltaQ(Y_bus,slackBusNumber,busData,Qspec,pqBuses,sBase,AngleData,V);
      deltaQ_compare = find_deltaQforPQBus(deltaQ,pqBuses);
      if sum(abs(deltaQ_compare)<epsilon) == length(deltaQ_compare)
          fprintf("deltaQ Converges\n")
          if sum(abs(deltaP(2:end))<epsilon) == length(deltaP)-1
              fprintf("deltaP also converges\n")
              k = 999;
          end
      else

          V = FDLF_Voltage(bPrime,bDoublePrime,deltaQ,deltaP,AngleData,V,slackBusNumber,pqBuses);
          Vold = [Vold V];
      end
      deltaP_old =[deltaP_old deltaP];
      deltaQ_old =[deltaQ_old deltaQ];
     fprintf("iteration :%d\n",k)

     k = k+1;

  end
  Angle = AngleData*(360/(2*pi));
  Pgen = P(1);
  Qgen = Q(1);
  Pgen_matrix = zeros(length(deltaP),1);
  Qgen_matrix = zeros(length(deltaQ),1);
  Pgen_matrix(1) = P(1);
  Qgen_matrix(1) = Q(1);
  for i = 1:length(pvBuses)
    if busData(pvBuses(i),4)>0
        Pgen = Pgen+ P(pvBuses(i)) + busData(pvBuses(i),4)/sBase;
        Pgen_matrix(pvBuses(i)) = P(pvBuses(i)) +busData(pvBuses(i),4)/sBase;
    else
        Pgen = Pgen+ P(pvBuses(i));
        Pgen_matrix(pvBuses(i)) = P(pvBuses(i));
    end

    if busData(pvBuses(i),5)>0
      Qgen = Qgen+ Q(pvBuses(i)) +busData(pvBuses(i),5)/sBase;
      Qgen_matrix(pvBuses(i)) = Q(pvBuses(i)) +busData(pvBuses(i),5)/sBase;
    else
      Qgen = Qgen+ Q(pvBuses(i));
      Qgen_matrix(pvBuses(i)) = Q(pvBuses(i));
    end
  end
  Pload = sum(busData(1:end,4))/sBase;
  Qload = sum(busData(1:end,5))/sBase;
  PowerLoss = Pgen-Pload;
  fprintf("Total Loss is %d MW",PowerLoss*sBase)
    
    
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

function [deltaP,P] = findDeltaP(Ybus,,Pspec,,Angle,V)
    P = zeros(length(Ybus),1);
    G = real(Ybus);
    B = imag(Ybus);
    for i = 1:length(Ybus)
        for k = 1:length(Ybus)
            P(i) = P(i) + V(k)*V(i)*((G(i,k)*cos(Angle(i)-Angle(k)))+(B(i,k)*sin(Angle(i)-Angle(k))));
        end
    end
    deltaP = Pspec-P;
end

function [deltaQ,Q] = findDeltaQ(Ybus,,,Qspec,,,Angle,V)
    Q = zeros(length(Ybus),1);
    G = real(Ybus);
    B = imag(Ybus);
    for i = 1:length(Ybus)
        for k = 1:length(Ybus)
            Q(i) = Q(i) + V(k)*V(i)*((G(i,k)*sin(Angle(i)-Angle(k)))-(B(i,k)*cos(Angle(i)-Angle(k))));
        end
    end

    deltaQ = Qspec-Q;

end

function [Angle] = FDLF_Angle(B1,B2,deltaP,Angle,V,slackBusNumber,pqBuses)

    V(slackBusNumber,:) = [];
    deltaP(slackBusNumber,:)= [];
    Angle(slackBusNumber,:)=[];
    P_process = deltaP./abs(V);
    deltaAngle = inv(B1)*P_process;
    Angle = Angle+deltaAngle;
    Angle = [0;Angle];

    for i=1:length(Angle)

        while Angle(i)>2*pi
            Angle(i) = Angle(i)-2*pi;
        end
        while Angle(i)<-2*pi
            Angle(i) = Angle(i)+2*pi;
        end
    end
end
function [V] = FDLF_Voltage(B1,B2,deltaQ,deltaP,Angle,V,slackBusNumber,pqBuses)
    V1 = V;
    V(slackBusNumber,:) = [];
    Vprocess = zeros(length(pqBuses),1);
    deltaQprocess = zeros(length(pqBuses),1);

    for i = 1:length(pqBuses)
        Vprocess(i) = V1(pqBuses(i));
        deltaQprocess(i) = deltaQ(pqBuses(i));
    end
    Qprocess = deltaQprocess./abs(Vprocess);
    deltaV = inv(B2)*Qprocess;
    Vprocess = Vprocess+deltaV;
    for i=1:length(Vprocess)
        V1(pqBuses(i)) = Vprocess(i);
    end

    V = V1;
end

function busData = finalVersionBusData(busData,slackBus,Angle,slackBusNumber)
    temp1 = busData(1:end,1:3);
    temp2 = busData(1:end,4:end);
    temp3 = slackBus(1:3);
    temp4 = slackBus(4:end);
    busData = [temp1,Angle,temp2];
    slackBus = [temp3,0,temp4];

    if slackBusNumber == 1
         busData = [slackBus;busData];
    else
        busData = [busData(1:slackBusNumber-1,1:end);slackBus;busData(slackBusNumber+1:end,1:end)];
    end
end



function Ybus = addDiagonal(Ybus1)
    Ybus = Ybus1;
    for i = 1:length(Ybus1)
        Ybus(i,i) = sum(Ybus1(i,1:end));
    end

end

function deltaQ_compare = find_deltaQforPQBus(deltaQ,pqBuses)
    temp = zeros(length(pqBuses),1);
    for i = 1:length(pqBuses)
        temp(i) = deltaQ(pqBuses(i));
    end
    deltaQ_compare = temp;

end

function [Pspec,Qspec]=  specValues(busData,slackBusNumber,sBase)

    Pspec = (busData(1:end,6)-busData(1:end,4))/sBase;
    Qspec = (busData(1:end,7)-busData(1:end,5))/sBase;
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

% function Bprime = createBprime(Y_bus,slackBusNumber)
%     B = -imag(Y_bus);
%     Bcopy = B;
%     
%     Bcopy(slackBusNumber,:) = [];
%     Bcopy(:,slackBusNumber) = [];
%     
%     Bprime = Bcopy * 1i;
% 
% end
% 
% 
% function BDoubleprime = createBDoublePrime(Y_bus, pqBuses,pvBuses,slackBus,bPrime,busNum,busData)
%       bPrimeCopy = bPrime;
%       pqBusNumber = length(pqBuses);
%       pvBusNumber = length(pvBuses);
%       slackBusNumber = findYBusEntries(slackBus,busNum,busData);
%       deleteCount = 0;
%       for i=1:pvBusNumber
%           numOfBus = findYBusEntries(pvBuses(i),busNum,busData);
%           
%           if numOfBus < slackBusNumber
%               bPrimeCopy(numOfBus-deleteCount,:) = [];
%               bPrimeCopy(:,numOfBus-deleteCount) = [];
%           else
%               bPrimeCopy(numOfBus-1-deleteCount,:) = [];
%               bPrimeCopy(:,numOfBus-1-deleteCount) = [];
%           end
%           deleteCount = deleteCount + 1;
%       end
%       
%       BDoubleprime = bPrimeCopy;
% end


function [B,bPrime,bDoublePrime] = findB(Y_bus,slackBusNum,pvBuses,pqBuses)
    B = imag(Y_bus);
    bPrime = -B;
    bPrime(:,slackBusNum) = [];
    bPrime(slackBusNum,:) = [];

    x = 1;
    y = 1;
    lenbDoublePrime = length(pqBuses);
    bDoublePrime = zeros(lenbDoublePrime,lenbDoublePrime);
    for i = 1:length(pqBuses)
        for k = 1:length(pqBuses)
            bDoublePrime(x,y) = -imag(Y_bus(pqBuses(k),pqBuses(i)));
            y=y+1;
        end
        x = x+1;
        y = 1;
    end
end
function [yBus] = writeOnBus(busData,branchData,busNum)
    yBus = zeros(busNum,busNum);
    % Find Yij entries of Y_bus matrix.
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

