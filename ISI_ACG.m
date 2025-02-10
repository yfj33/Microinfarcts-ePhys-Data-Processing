function [ISI,ACG] = ISI_ACG(ST)%SpikeTrain in s
%Output ISI ACG matrix for a given spike train, time 50ms/+-50ms with 0.5ms
%gap
ST=ST*1000; %convert to ms unit

ISI=diff(ST);
%ACG
 ACG=[];
for i=1:numel(ST)
    dis=ST-ST(i);
   
   ACG=[ACG dis];
end
ACG(ACG==0)=[];


end

