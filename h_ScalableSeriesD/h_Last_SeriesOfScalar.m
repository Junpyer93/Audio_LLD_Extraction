function lastValues = Last_SeriesOfScalar(auData, scalingRatio, elementNum, write_flag)

%%% File: Last_SeriesOfScalar.m
%%
%% --------------  ---------------
%%
%% The function of this subroutine is Series of Last 
%% of groups of samples. Size of the vector should equal 
%% elementNum, or 0 if Raw is present.
%% 
%% Definition:
%%
%% N = P*Q (?)
%% P -> scaleRatio, Q -> rescaleFactor 
%% LastValue = \Last_{i=1+(k-1)*N}^{kN} x_i
%% x -> audata
%% ignore when Weight present,  
%% (N3704)
%%
%%
%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)

global fid;

%%if nargin < 3, error('constr requires three input arguments'); end
%%if nargin < 5, weight_flag = 0; weight= []; end


%% Initialization
totalSampleNum = length(auData);
len = length(scalingRatio);

lastValues =  [];



if len == 1 & scalingRatio(len) == 1& write_flag == 1
	x = strcat('<Last/>');
   fprintf(fid, '%s\n ',x);
else
   lastlenData = 0;
   k = 1;
   len = length(scalingRatio);
   while k <= len
      lenData = lastlenData + scalingRatio(k)*elementNum(k);
      if totalSampleNum < lenData
      	N = floor((totalSampleNum - lastlenData)/scalingRatio(k));
         if N >= 0
       		for i = 1:N
         		sn_segment = auData(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
         		lastValues = [lastValues sn_segment(scalingRatio(k))];
            end
            %N1 = lastlenData + N*scalingRatio(k);
            %sn_segment = auData(N1+1:totalSampleNum);
            %LastValues = [LastValues, sn_segment(1)];
         end
      else
         for i = 1:elementNum(k)
            sn_segment = auData(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
 				lastValues = [lastValues sn_segment(scalingRatio(k))]
         end
         lastlenData = lenData;
      end
         if write_flag == 1
          	x = strcat('<Scaling ratio="',num2str(scalingRatio(k)),'" elementNum="',num2str(elementNum(k)),'"/>');
            fprintf(fid, '%s\n ',x);
         end
         k=k+1;
   end
   	x = strcat('<Last>', num2str(lastValues),'</Last>');
   	fprintf(fid,'%s\n',x);
end
 
    
       

   
