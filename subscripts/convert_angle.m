function out=convert_angle(in,rad)
%this function converts the rad values to deg and a scale ranging from 0 to
%360 
if nargin<2
ng=find(in<=0);
in(ng)=(in(ng))+2*pi;
out=in*180/pi;

elseif nargin==2 && strcmp(rad,'rad')
    
    ng=find(in<=0);
    in(ng)=(in(ng))+2*pi;
    out=in;
else
    disp('If you want to convert to rad put "rad" as second input, if you want deg, dont use any second input')
end 

end 