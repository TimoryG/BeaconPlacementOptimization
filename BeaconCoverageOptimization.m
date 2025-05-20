side_lengths = [36, 6, 16, 6, 36, 6, 16, 6];  % Side length here
lIncrement = 0.5; % Set the increment along the edges for placing beacons (length)
systematically_place_beacons(side_lengths, lIncrement);


function systematically_place_beacons(sides, lIncrement)
    % Function to systematically place beacons (circles) at all possible combinations
    % of three positions along polygon edges based on an incrementing variable.
    % sides: Array of side lengths for the polygon
    % increment: Length increment along each edge where beacons can be place
    setGlobalmC(0); % Set the maximum space of coverage intersection to 0
    setGlobalmCS(100000); % Set the minimum not covered space to 10000
    setGlobalmCT(0); % Set the maximum space of two beacon coverage to 0
    colorChoose = ['r', 'b', 'y']; % Set the color of the beacon ranges

    % Generate polygon vertices based on the specified side length
    [x, y] = create_custom_polygon(sides);
    num_edges = length(sides);

    % Find all combinations of edges 
    edge_combinations = nchoosek(1:num_edges, 3);  % All combinations of three edges
    beaconPoints = [0,0; 0,0; 0,0]; 

    % Loop through each combination of three edges
    for i = 1:size(edge_combinations, 1)
        edge_indices = edge_combinations(i, :);

        % Find all combinations of positions
        sideOne = 0:lIncrement/sides(edge_combinations(i, 1)):1;
        sideTwo = 0:lIncrement/sides(edge_combinations(i, 2)):1;
        sideThree = 0:lIncrement/sides(edge_combinations(i, 3)):1;
        pos_combinations = combvec(sideOne, sideTwo, sideThree); 

        % Loop through each combination of three positions along those edges
        for j = 1:size(pos_combinations, 2)
            position_points = pos_combinations(:, j); 

            %Uncomment these lines to see where the current position is
            %fprintf('New Position: '); 
            %fprintf('%d,%d,%d \n', position_points(1), position_points(2), position_points(3)), 

            % Create figure and plot the polygon
            figure;
            fill([x, x(1)], [y, y(1)], [0.5, 0.5, 0.5]);
            hold on;
            %Ensure the room stays the same number of pixels 
            axis([-50,100,-50,100]);
            set(gcf, "ToolBar", "none");
            set(gcf, "MenuBar", "none");
            set(gcf, "Position", [0, 0, 650, 600]); 

            % Place a beacon (circle) on each selected edge at the specified position
            for k = 1:3
                edge_index = edge_indices(k);
                position_point = position_points(k);
                % Get coordinates of the start and end points of the edge
                x1 = x(edge_index);
                y1 = y(edge_index);
                x2 = x(mod(edge_index, num_edges) + 1);
                y2 = y(mod(edge_index, num_edges) + 1);

                % Calculate the specified point along the edge
                point_x = x1 + position_point * (x2 - x1);
                point_y = y1 + position_point * (y2 - y1);
                beaconPoints(k,1)= point_x; 
                beaconPoints(k,2) = point_y; 

                % Draw a filled circle centered on the specified point
                theta = linspace(0, 2*pi, 100);
                circle_radius = 35;
                circle_x = point_x + circle_radius * cos(theta);
                circle_y = point_y + circle_radius * sin(theta);
                test = fill(circle_x, circle_y, colorChoose(k), 'EdgeColor', 'none', 'LineWidth', 1.5);
                test.FaceAlpha = 0.5; 

            end
            %Save the image and use imread to get the pixel breakdown
            exportgraphics(gcf, 'graph.png'); 
            img =imread('graph.png'); 

            %Calculate the amount of space covered by no beacons, two
            %beacons and three beacons
            coverage = findCoverageIntersection(img);
            coverageS = findNotCovered(img);
            coverageT = findCoverageTwo(img); 

            % Mark the center of the circle
            plot(beaconPoints(1,1), beaconPoints(1,2), 'ko', 'MarkerFaceColor', 'k');
            plot(beaconPoints(2,1), beaconPoints(2,2), 'ko', 'MarkerFaceColor', 'k');
            plot(beaconPoints(3,1), beaconPoints(3,2), 'ko', 'MarkerFaceColor', 'k');
            hold off

            % Check if beacons cover more of the room, if so update max
            % graph. 
            if coverageS < getGlobalmCS
                %Update best variables 
                setGlobalmC(coverage); 
                setGlobalmCS(coverageS);
                setGlobalmCT(coverageT);
                %Set the best position to the current one 
                bestPosition = [edge_indices(1) position_points(1) edge_indices(2) position_points(2) edge_indices(3) position_points(3)]; 
                setGlobalbP(bestPosition); 
                %Export the new best position to the maximalGraph png
                exportgraphics(gcf, 'maximalGraph.png'); 
            end 
            if coverageS == getGlobalmCS
            % Check if more area is covered by all three beacons
                if coverage > getGlobalmC
                    %Update best variables 
                    setGlobalmC(coverage); 
                    setGlobalmCS(coverageS); 
                    setGlobalmCT(coverageT);
                    %Set the best position to the current one 
                    bestPosition = [edge_indices(1) position_points(1) edge_indices(2) position_points(2) edge_indices(3) position_points(3)];
                    setGlobalbP(bestPosition); 
                    %Export the new best position to the maximalGraph png
                    exportgraphics(gcf, 'maximalGraph.png');
                end 
                if coverage == getGlobalmC
                % Check if more area is covered by two beacons 
                    if coverageT > getGlobalmCT
                        %Update best variables 
                        setGlobalmC(coverage); 
                        setGlobalmCS(coverageS);
                        setGlobalmCT(coverageT);
                        %Set the best position to the current one 
                        bestPosition = [edge_indices(1) position_points(1) edge_indices(2) position_points(2) edge_indices(3) position_points(3)];
                        setGlobalbP(bestPosition); 
                        %Export the new best position to the maximalGraph png
                        exportgraphics(gcf, 'maximalGraph.png');
                    end 
                    if coverageT == getGlobalmCT
                        %Add the position to the best position vector
                        rowVectorToInsert = [edge_indices(1) position_points(1) edge_indices(2) position_points(2) edge_indices(3) position_points(3)];
                        bestPositionNew = [bestPosition(1:end,:); rowVectorToInsert];
                        bestPosition= bestPositionNew; 
                        setGlobalbP(bestPosition); 
                    end 
                end
            end 
            close(gcf);
        end
    end
    %Reformat the best position matrix for output
    FinalResult = transpose(getGlobalbP);
    fprintf('Best Position=')
    fprintf([repmat(' %1.11f', 1, size(FinalResult, 1)) '\n'], FinalResult(:,:)) 
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

function [numOfPixels] = findCoverageIntersection(img)
    %Function to calculate how much of the room is covered by three beacons
    target_rbg = [175, 143, 79]; 
    mask = img(:,:,1) == target_rbg(1) & img(:,:,2) == target_rbg(2) & img(:,:,3) == target_rbg(3);
    numOfPixels = sum(mask(:));
end 

function [numOfPixels] = findCoverageTwo(img)
    %Function to calculate how much of the room is covered by two beacons
    target_rbg = [159,159,95; 95,31,159; 223,159,31];
    numOfPixels =0; 
    for i = 1:3
    mask = img(:,:,1) == target_rbg(i,1) & img(:,:,2) == target_rbg(i,2) & img(:,:,3) == target_rbg(i,3);
    numOfPixels = sum(mask(:));
    end 
end

function [numOfPixels] = findNotCovered(img)
    %Function to calculate how much of the room is covered by no beacons
    target_rbg = [128, 128, 128]; 
    mask = img(:,:,1) == target_rbg(1) & img(:,:,2) == target_rbg(2) & img(:,:,3) == target_rbg(3);
    numOfPixels = sum(mask(:));
end

function setGlobalmC(val)
    global maxCoverage
    maxCoverage = val; 
end

function r = getGlobalmC
    global maxCoverage
    r = maxCoverage; 
end

function setGlobalmCS(val)
    global maxCoverageS
    maxCoverageS = val; 
end

function r = getGlobalmCS
    global maxCoverageS
    r = maxCoverageS; 
end

function setGlobalmCT(val)
    global maxCoverageT
    maxCoverageT = val; 
end

function r = getGlobalmCT
    global maxCoverageT
    r = maxCoverageT; 
end

function setGlobalbP(val)
    global bestPosition
    bestPosition = val; 
end

function r = getGlobalbP
    global bestPosition
    r=bestPosition; 
end


