clc; clearvars; close all;

%% Initialization
S = randi([1 200], 3, 200);
%figure,
%scatter3(S(1,1:9),S(2,1:9),S(3,1:9),'*','r');

Scount = zeros(1,200);      % Number of people having a particular skill
J = zeros(3,200);           % J contains all the possible job vectors
Jcount = zeros(1,200);      % Number of people having a particular job
MarkJ = zeros(1,200);       % A flag indicator for jobs
A = randi([20 30], 3, 3);   % Transformation matrix

%% DATA READER MODULE
[ID,T1] = xlsread('Data_1.xlsx');
JobTitle = T1(2:length(T1),3);
Skills = T1(2:length(T1),2);
MapCnt=1;JobMapCnt=1;
Names = cell(1);
map = containers.Map();    % HashMap for skillset -> 'KeyType',kType,'ValueType',vType
JobMap = containers.Map();                                  % HashMap for Job

%% Algorithm for training
for k=1:length(ID)              % For each user profile
    Ui = zeros(1,1);            % Create skill set of particular user
    str = char(Skills(k));
    str1 = char('');
    for i=1:length(str)
        if str(i)==','
            if map.size == 0
                map(str1) = MapCnt;
                Names(MapCnt) = cellstr(str1);
                MapCnt=MapCnt+1;
            else
                tf = isKey(map,str1);
                if tf==0 % Map doesn't have it
                    map(str1) = MapCnt;
                    Names(MapCnt) = cellstr(str1);
                    MapCnt=MapCnt+1;
                end
            end
            if Ui(1)==0
                Ui(1) = map(str1);
            else
                Ui = [Ui map(str1)];
            end
            str1='';i=i+1;
        else
            str1=strcat(str1,str(i));
        end
    end
    tf = isKey(map,str1);
    if tf==0                    % Map doesn't have it
        map(str1) = MapCnt;
        Names(MapCnt) = cellstr(str1);
        MapCnt=MapCnt+1;
    end
    if Ui(1)==0
        Ui(1) = map(str1);
    else
        Ui = [Ui map(str1)];
    end
    
    str = char(JobTitle(k));
    if JobMap.size == 0
        JobMap(str) = JobMapCnt;
        JobMapCnt=JobMapCnt+1;
    else
        tf = isKey(JobMap,str);
        if tf==0                % Map doesn't have it
            JobMap(str) = JobMapCnt;
            JobMapCnt=JobMapCnt+1;
        end
    end
    if tf==0
        UJi = JobMapCnt-1;
    else
        UJi = JobMap(str);
    end
    % Update all the skill vectors in Ui
    meanS = zeros(3,1);
    totCount=0;
    for it=1:length(Ui)
        meanS = meanS + (Scount(Ui(it))+1)*S(:,Ui(it));
        totCount = totCount + Scount(Ui(it))+1;
    end
    meanS = meanS/totCount;         % Calculating the mean vector
    
    % Updating the skill vectors
    for it=1:length(Ui)
        %Cross weights taken into account
        S(:,Ui(it)) = ( (totCount - Scount(Ui(it)) - 1) * S(:,Ui(it)) + (Scount(Ui(it))+1)*meanS)/(totCount);		
    end
    %Building the skill matrix Si for user i
    Si = zeros(3,length(Ui));
    for it=1:length(Ui)
        Si(:,it) = S(:,Ui(it));
        Scount(Ui(it))=Scount(Ui(it))+1;			%Also updating the Skill count
    end
    Bi = A*Si;          %Computing Bi, the matrix that is approximation to possible jobs
    if MarkJ(UJi)==0	% Job is new to the program
        MarkJ(UJi)=1;
        J(:,UJi) = mean(Bi')';
        Jcount(UJi) = 1;
    else				%Job is known to the program
        % Finding nearest vector in Bi that matches J(Uji)
        min = Inf;
        minI = 1;
        for j=1:length(Bi)
            if min>norm(Bi(j)-J(UJi))
                min=norm(Bi(j)-J(UJi));
                minI = j;
            end
            %Adjust the job vector according to weighted mean
            J(:,UJi) = (Jcount(UJi)*J(:,UJi) + Bi(minI))/(Jcount(UJi)+1);
            Jcount(UJi)=Jcount(UJi)+1;
        end
    end
    % disp(k);
    % disp(meanS);
    % disp(Ui);
end
%% Plotting for visualization
figure;
sct1 = scatter3(S(1,1:MapCnt-1),S(2,1:MapCnt-1),S(3,1:MapCnt-1),'filled','b');    hold on;

%% Finding other related skills
[ID2,T2] = xlsread('D1N.xlsx');     % Reading test Dataset
JobTitle2 = T2(2:length(T2),3);
Skills2 = T2(2:length(T2),2);

for k=1:length(ID2) % For each query
    Ui = zeros(1,1);
    % Get the corresponding skill's value
    str = char(Skills2(k));
    str2 = '';
    for i=1:length(str)
        if str(i)==','
            if Ui(1)==0
                Ui(1) = cell2mat(values(map,cellstr(str2)));
            else
                val = cell2mat(values(map,cellstr(str2)));
                Ui = [Ui val];
            end
            str2='';i=i+1;
        else
            str2=strcat(str2,str(i));
        end
    end
    if Ui(1)==0
        Ui(1) = cell2mat(values(map,cellstr(str2)));
    else
        val = cell2mat(values(map,cellstr(str2)));
        Ui = [Ui val];
    end
    % Create a Skill matrix X
    X = Ui;
    Xnew = zeros(1, 2*size(X,2));	%For each vector in X, we find 2 nearest skill vectors
    for i=1:size(X,2)
        %find nearest 2 skills
        min1 = Inf;
        min1skill = 1;
        min2 = Inf;
        min2skill = 2;
        for j=1:length(Names)
            if abs(norm(S(:,X(:,i))-S(:,j)))<min1 && j~=X(:,i)
                min2 = min1;
                min2skill = min1skill;
                min1 = abs(norm(S(:,X(:,i))-S(:,j)));
                min1skill = j;
            else if abs(norm(S(:,X(:,i))-S(:,j)))<min2 && j~=X(:,i)
                    min2 = abs(norm(S(:,X(:,i))-S(:,j)));
                    min2skill = j;
                end
            end
            Xnew(1,2*i-1) = min1skill;
            Xnew(1,2*i) =  min2skill;
        end
    end
    %Sorting and displaying "unique" results
    sort(Xnew);
    fprintf('Suggested skills for user %d:\n', k);
    disp(Names(Xnew(1))); disp(Names(Xnew(2)));
    S_sug = zeros(3,200); cnt=2;
    %disp(S(:,Xnew(1)));
    S_sug(:,1) = S(:,Xnew(1));
    for i=2:size(Xnew, 2)
        if Xnew(i)~=Xnew(i-1)
            %disp(S(:,Xnew(i)));
            S_sug(:,cnt) = S(:,Xnew(i));    cnt=cnt+1;
        end
    end
end
% Plotting for visualization
sct2 = scatter3(S_sug(1,1:MapCnt-1),S_sug(2,1:MapCnt-1),S_sug(3,1:MapCnt-1),'*','r');
title('Skill Space');
legend([sct1, sct2], {'Skill Vectors', 'Suggested Skills'});