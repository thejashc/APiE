% Parameters
res = 50; % resolution for your spheres
rad = 0.25; % radius of spheres

% Make base sphere data
[xb, yb, zb] = sphere(res);

% Plot spheres
for i=1:5:length(t)
    
    % position of r
    surf(rad*xb + xA(i, 1), rad*yb + xA(i, 2), rad*zb + xA(i,3) , 'facecolor', 'r', 'edgealpha',0)
    hold on;
    surf(rad*xb + xB(i, 1), rad*yb + xB(i, 2), rad*zb + xB(i, 3), 'facecolor' ,'r', 'edgealpha',0)
    
    % shading
    light;
    lighting gouraud;
    
    axis square;
    xlim([-2 3])
    ylim([-2 3])
    zlim([-2 3])
    pause(0.01)
    
    hold off;
end

