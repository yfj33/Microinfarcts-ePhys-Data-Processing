function  plot_fill_circle(center,diameter,clr)
%plot a filled circle with diameter defined and color defined
%input: center of circle, diameter, and facecolor

%viscircles(center,diameter/2,'Edgecolor','none'); 
%fill circle
th=linspace(0,2*pi,100);
circ_x2=center(1)+diameter/2*cos(th);
circ_y2=center(2)+diameter/2*sin(th);
fill(circ_x2,circ_y2,clr,'EdgeColor','none','FaceAlpha',0.27);

end