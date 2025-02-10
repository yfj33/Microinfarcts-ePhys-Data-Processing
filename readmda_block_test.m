clc;clear;close all
time=5;
filt='filt.mda';
firings='firings.mda';
ts_curation='ts_curation.mat';
sr=30000;
%firings=readmda('firings.mda');
firings=readmda(firings);
PC=firings(1,:);
Event=firings(2,:);
Clust=firings(3,:);
ts_curation=load(ts_curation);
ts_curation=ts_curation.ts_curation;
ClustNum=numel(ts_curation);

if (strcmp(filt(end-4:end),'.csv')==1)
    A=textread(filt,'','delimiter',',');
    return;
end


F=fopen(filt,'rb');

try
code=fread(F,1,'int32');%code<0
catch
    error('Problem reading file: %s',filt);
end
if (code>0) 
    num_dims=code;
    code=-1;
else
    fread(F,1,'int32');%read dimension number
    num_dims=fread(F,1,'int32');    
end;

dim_type_str='int32';
if (num_dims<0)
    num_dims=-num_dims;
    dim_type_str='int64';
end;

S=zeros(1,num_dims);
for j=1:num_dims
    S(j)=fread(F,1,dim_type_str);
end;

%storeed data number in this time block
block=time*60*sr; 
% N=prod(S);

%N=S(1)*block;

% if num_dims == 1;
%   A = zeros(1,block);
% else
%   A=zeros(S(1),block);
% end
index=0;
indx=0;
Points_to_load=S(2);
pts=Points_to_load;

templates_new=zeros(S(1),100,ClustNum);%preallocate matrix to store template
Events_count=zeros(S(1),ClustNum);%preallocate matrix to store cluster count 
%while ~feof(F) %calculate templates with 5min as a block

for f=1: floor(pts/block)
    N=S(1)*block;
    A=zeros(S(1),block);
% else
%     N=S(1)*pts;%outside of last 5min
%     A=zeros(S(1),pts);

    index=index+time;
    indx=indx+1;
if (code==-1)
    M=zeros(1,N*2);
    M(:)=fread(F,N*2,'float');
    A(:)=M(1:2:prod(S)*2)+i*M(2:2:prod(S)*2);
elseif (code==-2)
    A(:)=fread(F,N,'uchar'); %5min of data
elseif (code==-3)
    A(:)=fread(F,N,'float');
elseif (code==-4)
    A(:)=fread(F,N,'int16');
elseif (code==-5)
    A(:)=fread(F,N,'int32');
elseif (code==-6)
    A(:)=fread(F,N,'uint16');
elseif (code==-7)
    A(:)=fread(F,N,'double');
elseif (code==-8)
    A(:)=fread(F,N,'uint32');
else
    error('Unsupported data type code: %d',code);
 disp(fprintf('Processing minute: %d',index))

end;

%[Template] = PCwaveform(PCevents,A(PCchannelnum),100);
%templates=[templates;Template];
interval_lo=(indx-1)*time*60*sr;
interval_up=(indx)*time*60*sr;

% event_of_interest=Event(Event>interval_lo & Event<interval_up);
% event_of_interest(1)=[];
% event_of_interest(end)=[];
% temp=cell(S(1),numel(event_of_interest));%store waveform for averaging
for n=1:ClustNum
 events=ts_curation{n,1};
 events_down=downsample(events,5); %down sample by 5 to get events number
 event_of_interest=events_down(events_down>interval_lo & events_down<interval_up);

 event_of_interest(1:3)=[];%get rid of first three
 event_of_interest(end-3:end)=[];%get rid of last three
 for t=1:numel(event_of_interest)
    for b=1:S(1)%each channel
    whpt= event_of_interest(t); %where event occurs
    start_point=[whpt-block*(indx-1)]-49;
    end_point=[whpt-block*(indx-1)]+50;
    templates_new(b,:,n)=A(b,start_point:end_point)+ templates_new(b,:,n); %need to substract previous session

    Events_count(b,n)= Events_count(b,n)+1;
    end
 end
end

% pts=pts-block;
end

Templates_New=zeros(S(1),100,ClustNum);
for j=1:S(1)
    for l=1:ClustNum
   Templates_New(j,:,l)=templates_new(j,:,l)./ Events_count(j,l);
    end
end

fclose(F);

 figure
for i=1:20
subplot(1,20,i);
plot(Templates_New(i,:,3));
ylim([-200 200]);
end
ylim([-100 100])