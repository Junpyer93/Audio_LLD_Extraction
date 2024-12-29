function SeriesOfScalarBegin(totalSampleNum)

%%% File: SeriesOfScalarBegin.m
%%
%% --------------  ---------------
%%
%% The function of this subroutine is to write series of scalar
%% such like raw, min,max,mean,variance,random,first,last into
%% XML file. 
%% 
%% <!-- ##################################################################### -->
%% <!-- Definition of abstract SeriesOfScalarType                             -->
%%  <!-- ##################################################################### -->
%%  <complexType name="SeriesOfScalarType" base="mpeg7:ScalableSeriesType"
%%               derivedBy="extension">
%%    <element name="Raw" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="Min" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="Max" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="Mean" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="Random" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="First" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="Last" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="Variance" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%    <element name="Weight" type="mpeg7:FloatVectorType" minOccurs="0"/>
%%  </complexType>
%%
%% (N3704)
%%
%%
%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)

%if nargin < 5, error('constr requires five input arguments'); end
%if nargin < 6, weight_flag = 0; weight= []; end
global fid;

fid = fopen('SeriesOfScalar.xml','w');
x='<!-- ##################################################################### -->';
fprintf(fid, '%s\n',x);
x= '<!-- Definition of abstract SeriesOfScalarType                            -->';
fprintf(fid, '%s\n',x);
x= '<!-- ##################################################################### -->';
fprintf(fid, '%s\n',x);
x='<complexType name="SeriesOfScalarType" base="mpeg7:ScalableSeriesType" derivedBy="extension">';
fprintf(fid, '%s\n',x);
x='<element name="Raw" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="Min" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="Max" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="Mean" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="Random" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="First" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="Last" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="Variance" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='<element name="Weight" type="mpeg7:FloatVectorType" minOccurs="0"/>';
fprintf(fid, '%s\n',x);
x='</complexType>';
fprintf(fid, '%s\n',x);

x = '<AudioSegment idref="101">';
fprintf(fid, '%s\n',x);
x = ' <!-- MediaTime, etc. -->';
fprintf(fid, '%s\n',x);
% <AudioPower>
x = '<HopSize timeunit="PT1N44000F">1</HopSize>';
fprintf(fid,'%s\n',x);
x = strcat('<Value totalSampleNum="', num2str(totalSampleNum), '">');
fprintf(fid,'%s\n',x);
