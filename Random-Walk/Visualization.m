%f_open = fopen('outputL40.txt', 'r');
f_open = fopen('output.txt', 'r');
dist_av=0;
N = 3;
T = 100;

% Initialize a matrix to hold the vector data
M = zeros(N, T);
%A_aver=zeros(1,T);

% Read the rest of the file, storing the data in the matrix
for i = 1:N
    for j = 1:T
        M(i, j) = fscanf(f_open, '%f', 1);
    end
end
% Close the file
fclose(f_open);

A=M';

% Create a figure
figure;
% Create an animation object
%ani = animation('Position',[0 0 560 420]);

% Loop through the coordinates
for i = 1:100
    % Get the current coordinate
    x = A(i,1);
    y = A(i,2);
    z = A(i,3);
    
    % Plot the point
    scatter3(x, y, z);
    hold on
    
    % Set the axis limits
    xlim([0 20]);
    ylim([0 20]);
    zlim([0 20]);
    
    % Set the view angle
    view(2);
    % Add the frame to the animation
    % ani.Frame(i) = getframe(gca);
    % Pause for a certain amount of time
    pause(1);
    
    % Clear the previous point
    clf;
end
%movie2avi(ani, 'animation.avi')

% Plot the vectors

% plot(t, M(1,:),'-r') %only one random walk
% hold on
%plot(t, A_aver, '-b') %avarage random walk
%legend({'Average distance'},'Location','northwest')
%title('Random walk in 3D Lattice with size of 40x40x40')
%xlabel('Time steps t') 
%ylabel('Distance values')
