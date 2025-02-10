function add_scalebar(time_scale,volt_scale,sr,lowestx,lowesty,ylim,xlim)
%add scalebar for waveform
%y: 1space is 1uV and x: 1space 1/sr s
Ylim=ylim(2)-ylim(1);
Xlim=xlim(2)-xlim(1);
xy_con=Ylim/Xlim;
con_fac=sr*10^-3;%ms
scalelength_x=time_scale*con_fac;
scalelength_y=0.11*scalelength_x;
scaleheight_x=10;
scaleheight_y=volt_scale;
scale_pos_x=[lowestx,lowesty,scalelength_x,scaleheight_x];
scale_pos_y=[lowestx,lowesty,scalelength_y,scaleheight_y];
rectangle('Position',scale_pos_x,'FaceColor','k','EdgeColor','k');
rectangle('Position',scale_pos_y,'FaceColor','k','EdgeColor','k');
%axis off
end