function [value,isterminal,direction] = myEventsFcn(t,y,p)
    % Locate the time when height passes through zero in a decreasing
    % direction
    % and stop integration.
    value = [y(1)-1; y(2)-1];
    isterminal = [1; 1];
    direction = [-1; -1];
end