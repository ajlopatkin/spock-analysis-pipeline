function bool = in_circle(center, radius, x, y)
%Check whether or not point (x,y) lies in the circle defined by inputs of
%center and radius
    bool = (x-center(:,1)).^2+(y-center(:,2)).^2<=radius.^2;
end