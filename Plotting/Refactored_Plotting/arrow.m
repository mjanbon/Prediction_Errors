function arrow(start_point, end_point, varargin)
    % Simple arrow drawing function if the arrow function is not available
    p = inputParser;
    addParameter(p, 'Length', 10);
    addParameter(p, 'TipAngle', 20);
    parse(p, varargin{:});
    
    % Draw the line
    line([start_point(1), end_point(1)], [start_point(2), end_point(2)], 'Color', 'k');
    
    % Calculate arrow head
    angle = atan2(end_point(2) - start_point(2), end_point(1) - start_point(1));
    tip_angle = p.Results.TipAngle * pi/180;
    length = p.Results.Length / 100;
    
    % Draw arrow head
    x1 = end_point(1) - length * cos(angle + tip_angle);
    y1 = end_point(2) - length * sin(angle + tip_angle);
    x2 = end_point(1) - length * cos(angle - tip_angle);
    y2 = end_point(2) - length * sin(angle - tip_angle);
    
    fill([end_point(1), x1, x2], [end_point(2), y1, y2], 'k');
end