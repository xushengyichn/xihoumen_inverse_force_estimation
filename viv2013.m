function [result]=viv2013(n,displayDates)
    if nargin <2
        displayDates = true;
    end

    switch n
        case 1
            % vibac2_name=["2013-01-14 00-vibac2.txt";"2013-01-14 01-vibac2.txt"];
            % vibac3_name=["2013-01-14 00-VIBac3.txt";"2013-01-14 01-VIBac3.txt"];
            % vibac4_name=["2013-01-14 00-VIBac4.txt";"2013-01-14 01-VIBac4.txt"];
            startDate = datetime(2013,1,14,0,0,0);
            endDate = datetime(2013,1,14,2,0,0);
        case 2
            startDate = datetime(2013,1,16,3,0,0);
            endDate = datetime(2013,1,16,5,0,0);
        case 3
            startDate = datetime(2013,1,22,2,0,0);
            endDate = datetime(2013,1,22,4,0,0);
        case 4
            startDate = datetime(2013,2,6,0,0,0);
            endDate = datetime(2013,2,6,3,0,0);
        case 5
            startDate = datetime(2013,3,6,18,0,0);
            endDate = datetime(2013,3,6,20,0,0);
        case 6
            startDate = datetime(2013,4,3,15,0,0);
            endDate = datetime(2013,4,3,17,0,0);
        case 7
            startDate = datetime(2013,4,8,15,0,0);
            endDate = datetime(2013,4,8,16,0,0);
        case 8
            startDate = datetime(2013,4,10,19,0,0);
            endDate = datetime(2013,4,10,22,0,0);
        case 9
            startDate = datetime(2013,5,19,11,0,0);
            endDate = datetime(2013,5,19,13,0,0);
        case 10
            startDate = datetime(2013,5,21,18,0,0);
            endDate = datetime(2013,5,21,20,0,0);
        case 11
            startDate = datetime(2013,6,10,4,0,0);
            endDate = datetime(2013,6,10,7,0,0);
        case 12
            startDate = datetime(2013,7,9,17,0,0);
            endDate = datetime(2013,7,9,18,0,0);
        case 13
            startDate = datetime(2013,7,11,18,0,0);
            endDate = datetime(2013,7,11,19,0,0);
        case 14
            startDate = datetime(2013,7,16,15,0,0);
            endDate = datetime(2013,7,16,17,0,0);
        case 15
            startDate = datetime(2013,7,21,16,0,0);
            endDate = datetime(2013,7,21,20,0,0);
        case 16
            startDate = datetime(2013,8,1,12,0,0);
            endDate = datetime(2013,8,1,16,0,0);
        case 17
            startDate = datetime(2013,8,2,22,0,0);
            endDate = datetime(2013,8,3,0,0,0);
        case 18
            startDate = datetime(2013,8,10,22,0,0);
            endDate = datetime(2013,8,10,23,0,0);
        case 19
            startDate = datetime(2013,8,16,18,0,0);
            endDate = datetime(2013,8,16,19,0,0);
        case 20
            startDate = datetime(2013,8,24,15,0,0);
            endDate = datetime(2013,8,24,21,0,0);
        case 21
            startDate = datetime(2013,8,28,14,0,0);
            endDate = datetime(2013,8,28,17,0,0);
        case 22
            startDate = datetime(2013,8,29,14,0,0);
            endDate = datetime(2013,8,29,15,0,0);
            
    end

    timeDifference = endDate-startDate ;
    numberOfHours = hours(timeDifference);
    startDate.Format = 'yyyy-MM-dd HH';
    endDate.Format = 'yyyy-MM-dd HH';
    for i=1:numberOfHours
        vibac2_name(i)=char(startDate+hours(i-1))+"-vibac2.txt";
        vibac3_name(i)=char(startDate+hours(i-1))+"-VIBac3.txt";
        vibac4_name(i)=char(startDate+hours(i-1))+"-VIBac4.txt";
    end


    if displayDates
        disp('Start Date and Time:');
        disp(startDate);
        disp('End Date and Time:');
        disp(endDate);
    end

    result.vibac2_name=vibac2_name;
    result.vibac3_name=vibac3_name;
    result.vibac4_name=vibac4_name;
    result.startDate=startDate;
    result.endDate=endDate;
end