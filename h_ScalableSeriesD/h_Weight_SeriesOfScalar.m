function weightValues = Weight_SeriesOfScalar(scalingRatio, elementNum, weight, write_flag)

%%% File: Weight_SeriesOfScalar.m
%%
%% --------------  ---------------
%%
%% The function of this subroutine is Series of Weight 
%% of groups of samples. Size of the vector should equal 
%% elementNum, or 0 if Raw is present.
%% 
%% Definition:
%%
%% N = P*Q (?)
%% P -> scaleRatio, Q -> rescaleFactor 
%% WeightValue = \sum_{i=1+(k-1)*N}^{kN} x_i
%% x -> audata
%% If Weight present, ignore samples with zero weight. 
%% If all have zero weight, set to zero by convention. 
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
totalSampleNum = length(weight);
len = length(scalingRatio);

weightValues =  [];



if len == 1 & scalingRatio(len) == 1 & write_flag == 1
    x = strcat('<Weight/>');
    fprintf(fid, '%s\n ',x);
else
   if sum(weight) == 0
    	weightValues = [weightValues 0]
   else
      lastlenData = 0;
      k = 1;
      len = length(scalingRatio(k));
      while k <= len
      	lenData = lastlenData + scalingRatio(k)*elementNum(k);
         if totalSampleNum < lenData
        		N = floor((totalSampleNum - lastlenData)/scalingRatio(k));
           	if N >= 0
       			for i = 1:N
                  weight_segment = weight(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
                  weightValues = [weightValues mean(weight_segment)];
            	end
              	N1 = lastlenData + N*scalingRatio(k);
               weight_segment = weight(N1+1:totalSampleNum);
               weightValues = [weightValues mean(weight_segment)];
           	end
        	else
           	for i = 1:elementNum(k)
                weight_segment =weight(lastlenData+1+(i-1)*scalingRatio(k):lastlenData+i*scalingRatio(k));
                weightValues = [weightValues mean(weight_segment)];
            end
            lastlenData = lenData;
         end
            if write_flag == 1
          		x = strcat('<Scaling ratio="',num2str(scalingRatio(k)),'" elementNum="',num2str(elementNum(k)),'"/>');
               fprintf(fid, '%s\n ',x);
            end
          	k=k+1;
      end
     end;
       x = strcat('<Weight>', num2str(weightValues),'</Weight>');
       fprintf(fid,'%s\n',x);
  end
    
 
 
 
    
       

   
