function g=plot_err_patch_v2(x,y,e,col1,col2,style)

x=x(:)'; %ms: inverts to column vector
y=y(:)';
e=e(:)';

ye1=y+e;
ye2=y-e; ye2=ye2(end:-1:1);
ye=[ye1,ye2];
xe=[x,x(end:-1:1)];

hold on;
h=patch(xe,ye,col2,'linestyle','none', 'FaceAlpha', 0.5);
if(nargin==6)
    g=plot(x,y,'color',col1,'linewidth',2,'linestyle',style);
else
    g=plot(x,y,'color',col1,'linewidth',2);
end