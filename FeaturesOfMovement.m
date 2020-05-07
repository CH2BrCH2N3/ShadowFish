clear all;

filestoprocess = dir('*.csv'); % stores details on all .csv files in the current folder

for f = 1 : length(filestoprocess) % loops for each .csv file in the folder
     
    clearvars -except filestoprocess f distTravelArray velocityArray timeMovingArray fileNames activeTimeArray actTailbeatFreqArray lowTailbeatFreq medTailbeatFreq highTailbeatFreq avgTBangleArray lowPeakArray medPeakArray highPeakArray
    
    fileID = fopen(filestoprocess(f).name); % opens current file as fileID
    
    disp(filestoprocess(f).name) % displays the name of the file being processed

    C = textscan(fileID,'%d %d %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d','Delimiter',','); % stores values from current file in cell array

    % stores spine angles in an array and converts to degrees
    spine(1,:)=C{3};
    spine(2,:)=C{4};
    spine(3,:)=C{5};
    spine(4,:)=C{6};
    spine(5,:)=C{7};
    spine = double(spine)*(180/pi);
    
    % stores X coordinates in an array and converts to cm
    nodeXCoord(1,:)=C{8};
    nodeXCoord(2,:)=C{10};
    nodeXCoord(3,:)=C{12};
    nodeXCoord(4,:)=C{14};
    nodeXCoord(5,:)=C{16};
    nodeXCoord(6,:)=C{18};
    nodeXCoord(7,:)=C{20};
    nodeXCoord = double(nodeXCoord)/1280*13.4;

    % stores Y coordinates in an array and converts to cm
    nodeYCoord(1,:)=C{9};
    nodeYCoord(2,:)=C{11};
    nodeYCoord(3,:)=C{13};
    nodeYCoord(4,:)=C{15};
    nodeYCoord(5,:)=C{17};
    nodeYCoord(6,:)=C{19};
    nodeYCoord(7,:)=C{21};
    nodeYCoord = double(nodeYCoord)/720*7.2;
    nodeYCoord = 7.2 - nodeYCoord; % corrects direction of y axis
    
    %% Calculates the distance travelled per second in cm using node 4
    
    n = 1;
        for i = 101:100:length(nodeXCoord)
            distSec(n) = sqrt(double((nodeXCoord(4,i)-nodeXCoord(4,i-100))^2+(nodeYCoord(4,i)-nodeYCoord(4,i-100))^2));
            n = n + 1;
        end
        
    distTravel = sum(distSec);
    velocity = distTravel/(length(distSec));
 
     %% Sums together angles for each frame
    
    for i = 1:length(spine)
        headtoTail(i) = 0;
        for j = 1:5
        if abs(spine(j,i)) > 5 
            headtoTail(i) = headtoTail(i) + abs(spine(j,i));
        else
        end
        end
    end 
    
    %% Plots the sum of angles over time
    
    smoothedTail = smooth(headtoTail,5); % running average over 5 frames
    
    figure;
    axisTime = [0.01:0.01:(length(spine)/100)];
    
    [pos_pks,pos_locs] = findpeaks(smoothedTail,'MinPeakDistance',5,'MinPeakProminence',5); % finds the peaks
    
    noSpeed = 0;
    lowSpeed = 0;
    medSpeed = 0;
    highSpeed = 0;
    
        for i = 1:length(distSec) % loops for each second
            
        j= i*100+1; % j is the counter for the final frame in a second
        k =j-100; % k is the counter for the first frame in a second
       
            if distSec(i) < 0.5 % counts less than 5mm in 1 second as no movement
                
               movement(i) = 0; % array for which seconds have movement and which don't
               
               noSpeed = noSpeed + 1; % counter for how many seconds the fish is stationary
               statPeaks = find(pos_locs>k & pos_locs<j); % finds the peaks in this second
               
               pos_pks(statPeaks)=NaN; % replaces the peaks with NaN in the main peak array
               pos_locs(statPeaks)=NaN;
               
        elseif distSec(i) < 2 % counts less than 2cm in 1 second as low speed
               movement(i) = 1; 
               lowSpeed = lowSpeed + 1; %low speed for this second
               
               lowPeaks = find(pos_locs>k & pos_locs<j);
               lowPeakAng = pos_pks(lowPeaks);
               
               if lowSpeed == 1
                    lowCount = vertcat(lowPeaks);
                    lowPeakAngles = vertcat(lowPeakAng);
               else
                   lowCount = vertcat(lowCount,lowPeaks);
                   lowPeakAngles = vertcat(lowPeakAngles,lowPeakAng);
               end
               
        elseif distSec(i) <= 4
               movement(i) = 1;
               medSpeed = medSpeed + 1; % medium speed
               
               medPeaks = find(pos_locs>k & pos_locs<j);
               medPeakAng = pos_pks(medPeaks);
               
               if medSpeed == 1
                    medCount = vertcat(medPeaks);
                    medPeakAngles = vertcat(medPeakAng);
               else
                   medCount = vertcat(medCount,medPeaks);
                   medPeakAngles = vertcat(medPeakAngles, medPeakAng);
               end
               
        else
               movement(i) = 1;
               highSpeed = highSpeed + 1; % high speed
               
               highPeaks = find(pos_locs>k & pos_locs<j);
               highPeakAng = pos_pks(highPeaks);
               
                if highSpeed == 1
                   highCount = vertcat(highPeaks);
                   highPeakAngles = vertcat(highPeakAng);
               else
                   highCount = vertcat(highCount,highPeaks);
                   highPeakAngles = vertcat(highPeakAngles, highPeakAng);
               end
        end
        end

        if lowSpeed == 0
            avgLowPeak = NaN;
        else
            avgLowPeak = mean(lowPeakAngles);
        end
        
        if medSpeed == 0
            avgMedPeak = NaN;
        else
            avgMedPeak = mean(medPeakAngles);
        end
  
        if highSpeed == 0
            avgHighPeak = NaN;
        else
            avgHighPeak = mean(highPeakAngles);
        end
        
        if lowSpeed>0
        lowTBfreq = length(lowCount)/lowSpeed;
        else
            lowTBfreq = NaN;
        end
        
        if medSpeed > 0
        medTBfreq = length(medCount)/medSpeed;
        else
        medTBfreq = NaN;
        end
        
        if highSpeed>0
        highTBfreq = length(highCount)/highSpeed;
        else
        highTBfreq = NaN;
        end
        
        % creates new array for peaks removing NaN from stationary seconds
        n = 1;
        for i=1:length(pos_locs)
            if isnan(pos_locs(i))
            else
                nanPos_locs(n) = pos_locs(i);
                nanPos_pks(n)= pos_pks(i);
                n = n+1;
            end
        end
    
        avgTailBendAngle = mean(nanPos_pks); %calculates avg TB angle only when the fish is moving
        
    plot(axisTime, smoothedTail,'-k');
    hold on
    plot(axisTime(nanPos_locs), nanPos_pks, 'or');
    hold off
    xlabel('Time (s)');
    ylabel('Tail bend amplitude');
    
    tailbeatFreq = length(nanPos_locs)/((length(spine)/100));
    
    timeMoving = ((sum(movement(:) == 1))/length(movement))*100;
    actTailbeatFreq = tailbeatFreq*(100/timeMoving);
    velocity = velocity * (timeMoving/100); % calculates velocity for when fish is actively moving
    
    activeTime = zeros(1,1000);
    
    m = 1;
    n = 1;
    for i = 1:length(movement)
        if movement(i) == 0 
            if i == length(movement)
            elseif movement(i+1) == 1
                m = m + 1; % if there is movement in the next 1s period we update the restTime counter
            else
            end
        elseif movement(i) == 1
            activeTime(n) = activeTime(n)+1; % with movement 1s is added to the active time duration
            if i == length(movement)
            elseif movement(i+1) == 0
                n = n + 1;  % if there is no movement in the next 1s period we update the activeTime counter 
            else
            end
        end
    end
    
    if activeTime(n) == 0
        n = n-1;
    else
    end
    
    activeTimeAvg = sum(activeTime)/n;
    
    %% Storing calculated values in arrays
    
    fileNames(f,1) = string(filestoprocess(f).name);
    distTravelArray(f,1) = distTravel;
    velocityArray(f,1) =  velocity;
    timeMovingArray(f,1) = timeMoving;
    activeTimeArray(f,1) = activeTimeAvg;
    actTailbeatFreqArray(f,1) = actTailbeatFreq;
    lowTailbeatFreq(f,1) = lowTBfreq;
    medTailbeatFreq(f,1) = medTBfreq;
    highTailbeatFreq(f,1) = highTBfreq;
    avgTBangleArray(f,1) = avgTailBendAngle;
    lowPeakArray(f,1) = avgLowPeak;
    medPeakArray(f,1) = avgMedPeak;
    highPeakArray(f,1) = avgHighPeak;

end

LowTBfreq = lowTailbeatFreq(~isnan(lowTailbeatFreq));
MedTBfreq = medTailbeatFreq(~isnan(medTailbeatFreq));
HighTBfreq = highTailbeatFreq(~isnan(highTailbeatFreq));
highSpeedPeak = highPeakArray(~isnan(highPeakArray));
medSpeedPeak = medPeakArray(~isnan(medPeakArray));
lowSpeedPeak = lowPeakArray(~isnan(lowPeakArray));

%% Outputting data in a .csv file with averages

file = 'FishAnalyserData.csv'; % creates .csv file to save all FishAnalyser data
outputHeaders = ["File", "Distance Travelled(cm)", "Mean active velocity(cm s^-1)", "Mean active tail beat frequency" "Time spent moving(%)", "Mean active time duration(s)", "Low speed tail beat frequency", "Medium speed tail beat frequency", "High speed tail beat frequency", "Mean tail bend amplitude", "Mean low speed tail bend amplitude", "Mean medium speed tail bend amplitude", "Mean high speed tail bend amplitude"];
outputData = cell(length(filestoprocess),13);
outputData = [fileNames, distTravelArray, velocityArray, actTailbeatFreqArray, timeMovingArray, activeTimeArray, lowTailbeatFreq, medTailbeatFreq, highTailbeatFreq, avgTBangleArray, lowPeakArray, medPeakArray, highPeakArray];
averages = ["Mean", mean(distTravelArray), mean(velocityArray), mean(actTailbeatFreqArray), mean(timeMovingArray), mean(activeTimeArray), mean(LowTBfreq), mean(MedTBfreq), mean(HighTBfreq), mean(avgTBangleArray), mean(lowSpeedPeak), mean(medSpeedPeak), mean(highSpeedPeak)];
standardDeviations = ["Standard Deviation", std(distTravelArray), std(velocityArray), std(actTailbeatFreqArray), std(timeMovingArray), std(activeTimeArray), std(LowTBfreq), std(MedTBfreq), std(HighTBfreq), std(avgTBangleArray), std(lowSpeedPeak), std(medSpeedPeak), std(highSpeedPeak)];
    
fID = fopen('FishAnalyserData.csv','wt');
fprintf(fID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',outputHeaders(:));
for i = 1:(length(filestoprocess))
    fprintf(fID,'%s, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n',outputData(i,1:13));
end
fprintf(fID,'\n%s, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n', averages(:));
fprintf(fID,'\n%s, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f', standardDeviations(:));
fclose(fID);

fclose(fileID);
    
           