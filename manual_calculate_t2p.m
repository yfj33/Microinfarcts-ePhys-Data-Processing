function [amp,t2p_manual] = manual_calculate_t2p(waveform,SR)
%Calculate t2p value manually
WF=waveform.filtWaveform;
t2p_manual=zeros(1,numel(WF));
amp=zeros(1,numel(WF));
for i=1:numel(WF)
wf=WF{i};    
[tval,tindx]=max(wf);
[pval,pindx]=min(wf);
t2p_manual(1,i)=abs(pindx-tindx)/SR*1000; %in ms
amp(1,i)=abs(pval-tval);
end

end