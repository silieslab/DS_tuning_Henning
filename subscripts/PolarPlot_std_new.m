% Next you need to add error bars. It is a matter of drawing small lines in 
% radial direction around each average value. Let us say you want to draw an
% error bar at the 3rd data point. For this error bar (a radial line), i need
% at least two points which are the error limits. These points would be 
% {average(3)-error(3)} and {average(3)+error(3)}. To plot this error bar 
% (radial line), i need to keep the angle constant. That is my last step
% 
% 3) polar([angle(3) angle(3)],[average(3)-error(3), average(3)+error(3)]).
% Repeat this for all data points.
% 
% Hope this makes it clear. You can modify the code easily. 



for ii=1:length(theta)-1
    
    e=stuningPdir(ii);
    m=mtuningPdir(ii);
    angle2=theta(ii);
    P=polar([angle2 angle2], [m-e, m+e]);
    P.Color=Color(i,:);
    
end 