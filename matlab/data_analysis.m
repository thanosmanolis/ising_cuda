%% SEQ & CUDA implementation

for i = 1:4
    if i == 1
        file = '../ising_cuda/results/v1.txt';
        line = '-';
        color = 'r';
    elseif i == 2
        file = '../ising_cuda/results/v2.txt';
        line = '-';
        color = 'b';
    elseif i == 3
        file = '../ising_cuda/results/v3.txt';
        line = '-';
        color = 'g';
    else
        file = '../ising_cuda/results/seq.txt';
        line = '-';
        color = [0.8500, 0.3250, 0.0980];
    end


    A_temp = load(file);
    A = sortrows(A_temp);
    A_25 = A(A(:,2)==25, :);
    A_50 = A(A(:,2)==50, :);
    A_100 = A(A(:,2)==100, :);

    x_25 = A_25(:,1);
    y_25 = A_25(:,2);
    z_25 = A_25(:,3);

    x_50 = A_50(:,1);
    y_50 = A_50(:,2);
    z_50 = A_50(:,3);

    x_100 = A_100(:,1);
    y_100 = A_100(:,2);
    z_100 = A_100(:,3);

    plot3(x_25, y_25, z_25, line, 'color', color);
    hold on
    plot3(x_50, y_50, z_50, line, 'color', color);
    hold on
    plot3(x_100, y_100, z_100, line, 'color', color);
end

%% General plot stuff
title('Athanasios Manolis 8856 - Ising Model');
xlabel('n - Size of matrix (nxn)');
xticks([0 2000 4000 6000 8000 10000])
ylabel('k - Iterations');
yticks([25 50 100])
zlabel('time (sec)');
grid on;
axis square;

% TextBoxes
annotation('textbox', [.7 .4 .4 .3], 'String', 'Blue: V1', 'FitBoxToText', 'on', 'Color', 'Blue');
annotation('textbox', [.7 .5 .4 .3], 'String', 'Red: V2', 'FitBoxToText', 'on', 'Color', 'Red');
annotation('textbox', [.7 .5 .4 .3], 'String', 'Green: V3', 'FitBoxToText', 'on', 'Color', 'Green');
annotation('textbox', [.7 .5 .4 .3], 'String', 'Orange: V0', 'FitBoxToText', 'on', 'Color', [0.8500, 0.3250, 0.0980]);
