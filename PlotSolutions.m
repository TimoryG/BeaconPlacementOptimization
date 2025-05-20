side_lengths = [60,60,60,60];  % Side length here
edge_positions = [1, 0.666; 2, 1; 3, 0.666]; %Beacon positions 
plot_custom_polygon_with_specified_circles(side_lengths, edge_positions)

function plot_custom_polygon_with_specified_circles(sides, edge_positions)
    % Function to create and plot a polygon with circles centered at specific points along edges
    % sides: Array of side lengths for the polygon
    % edge_positions: Array of tuples [edge_index, position_fraction],
    % where each tuple specifies the edge and position along that edge (0 to 1) for each circle
    
    % Generate polygon vertices based on the specified side lengths
    [x, y] = create_custom_polygon(sides);
    colorChoose = ['r', 'b', 'y'];

    % Plot the polygon
    fill([x, x(1)], [y, y(1)], [0.5, 0.5, 0.5]);
    hold on;
    %Ensure the room stays the same number of pixels 
    axis([-50,100,-50,100]);
    set(gcf, "ToolBar", "none");
    set(gcf, "MenuBar", "none");
    set(gcf, "Position", [0, 0, 650, 600]);

    % Loop through each specified circle position
    for i = 1:size(edge_positions, 1)
        edge_index = edge_positions(i, 1);
        position_fraction = edge_positions(i, 2);

        % Get coordinates of the start and end points of the edge
        x1 = x(edge_index);
        y1 = y(edge_index);
        x2 = x(mod(edge_index, length(sides)) + 1);
        y2 = y(mod(edge_index, length(sides)) + 1);

        % Calculate the specified point along the edge
        point_x = x1 + position_fraction * (x2 - x1);
        point_y = y1 + position_fraction * (y2 - y1);

        % Draw the circle centered on the specified point
        theta = linspace(0, 2*pi, 100);
        circle_radius = 35;
        circle_x = point_x + circle_radius * cos(theta);
        circle_y = point_y + circle_radius * sin(theta);
        test = fill(circle_x, circle_y, colorChoose(i), 'EdgeColor', 'none', 'LineWidth', 1.5);
        test.FaceAlpha = 0.5;

        % Mark the center of the circle
        plot(point_x, point_y, 'ko', 'MarkerFaceColor', 'k');
    end
    %Export the graphic to a png
    exportgraphics(gcf, 'solution.png');
    hold off;

end

 

function [x, y] = create_custom_polygon(sides)
    % Function to create polygon vertices based on specified side lengths
    num_sides = length(sides);
    angles = linspace(0, 2 * pi, num_sides + 1);
    x = zeros(1, num_sides + 1);
    y = zeros(1, num_sides + 1);

    % Calculate vertices based on each side length and angle
    for i = 1:num_sides
        x(i + 1) = x(i) + sides(i) * cos(angles(i));
        y(i + 1) = y(i) + sides(i) * sin(angles(i));
    end

    % Remove the last redundant point to keep the polygon closed
    x = x(1:end-1);
    y = y(1:end-1);

end
