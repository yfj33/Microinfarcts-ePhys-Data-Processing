function [Template] = PCwaveform(events,filt,range)
%generate primary channel waveform matrix of given units
%Input: events:firing events matrix(1*x in dps) of PC channel, filt:filtwaveform matrix(1*x in dps) of primary channel,
%and template range(number of dps)
%Output: template matrix for given time duration of data, size of
%events*100pts

%template size is chnum*range (100 dps)
Template=zeros(numel(event1),range);
%grab the points at the range
for pts=1:numel(events)
    event_position=events(pts);
    start_point=event_position-(range/2-1);
    end_point=event_position+range/2;
    Template(pts,:)=filt(start_point:end_point); %Template with size of events*100 dps
end



