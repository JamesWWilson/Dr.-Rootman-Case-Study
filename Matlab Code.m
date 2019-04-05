
%% Part 1: Graph average responses of participants for eye 1 vs. expert. 
% Create table out of excel sheet of modified participants
rootman = readtable('rootman 1 winter 19 copy.xlsx');
rootman_coord = rootman(:,7:54);
rootman_coord2 = table;

[r,c] = size(rootman_coord);
for i=1:c
  rootman_coord2.(i) = str2double(rootman_coord{:,i});
end

% change the table to a matrix
rootman_mat = table2array(rootman_coord2);
rootman_mat(isnan(rootman_mat))=0;
matrices = zeros(1000,1000,r);

% for all of the participants create a matrix
for i = 1:r
    for j = 1:(c-1)
        if mod(j , 2) == 1
           x = rootman_mat(i,j);
           j = j+1;
           y = rootman_mat(i,j);
           if x ~= 0 && y~= 0
              matrices(x,y,i)= 1;
           end       
        end   
    end
end

matrices_summed = (sum(matrices,3))/r;

% Averaged - plot values to see what this looks like
I = mat2gray(transpose(matrices_summed));
imshow(I)
title('Averaged Clicks for all Participants'); 

% Try to fit a curve to the average response of participants
[row, col] = find(matrices_summed);
[p, S] = polyfit(row,col,2);
x1 = linspace(200,600,695);
[y1, delta] = polyval(p, x1, S);

figure
plot(row,col,'b.')
hold on
plot(x1,y1, 'r')
plot(x1,y1+2*delta,'m--',x1,y1-2*delta,'m--')
xlim([0 800])
ylim([0 800])
rectangle('position',[200 100 10 600])
title('Quadratic Fit of Data with 95% Prediction Interval')
legend('Data','Quadratic Fit','95% Prediction Interval')
caption = sprintf('y = %f * x^2 + %f * x + %f', p(1), p(2), p(3));
text(50, 750, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');

hold off

%This is not good enough, need a way to get rid of all of the outliers.

% Try going across eye by intervals of 10

%calculate the sd of each interval
std_colmns = zeros(44,1);
j = 1;
for i = [181 191 201 211 221 231 241 251 261 271 281 291 301 311 321 331 341 351 361 371 381 391 401 411 421 431 441 451 461 471 481 491 501 511 521 531 541 551 561 571 581 591 601 611]
    [xx, yy] = find(matrices_summed(i:(i+9),:));
    std_colmns(j) = std(yy);
    j=j+1;
end

std_colmns2 = std_colmns(~isnan(std_colmns));
mean_std = mean(std_colmns2);

%keep points which are 1 standard deviation from MEAN
matrices_summed_new = zeros(1000,1000);

med_yy = 0;
minusSD_yy = 0;
plusSD_yy = 0;

for i = [181 191 201 211 221 231 241 251 261 271 281 291 301 311 321 331 341 351 361 371 381 391 401 411 421 431 441 451 461 471 481 491 501 511 521 531 541 551 561 571 581 591 601 611]
    [xx, yy, z] = find(matrices_summed(i:(i+9),:));
    med_yy = mean(yy);
    minusSD_yy = med_yy - mean_std;
    plusSD_yy = med_yy + mean_std;
    for j = 1:length(yy)
        if (yy(j) < plusSD_yy) && (yy(j) > minusSD_yy)
            matrices_summed_new((xx(j) + i - 1), yy(j)) = z(j);
        end
    end    
end

%keep points which are 1 standard deviation from MEDIAN
matrices_summed_new_median = zeros(1000,1000);
med_yy = 0;
minusSD_yy = 0;
plusSD_yy = 0;

for i = [181 191 201 211 221 231 241 251 261 271 281 291 301 311 321 331 341 351 361 371 381 391 401 411 421 431 441 451 461 471 481 491 501 511 521 531 541 551 561 571 581 591 601 611]
    [xx, yy, z] = find(matrices_summed(i:(i+9),:));
    med_yy = median(yy);
    minusSD_yy = med_yy - mean_std;
    plusSD_yy = med_yy + mean_std;
    for j = 1:length(yy)
        if (yy(j) < plusSD_yy) && (yy(j) > minusSD_yy)
            matrices_summed_new_median((xx(j) + i - 1), yy(j)) = z(j);
        end
    end    
end

[row1, col1] = find(matrices_summed_new);
[row1med, col1med] = find(matrices_summed_new_median);

x1 = linspace(200,600,800);

[p, S] = polyfit(row1,col1,2);
[y1, delta] = polyval(p, x1, S);

[p_med, S_med] = polyfit(row1med,col1med,2);
[y1_med, delta_med] = polyval(p_med, x1, S_med);

% graph mean vs median
figure
plot(row1,col1,'b.')
hold on
plot(x1,y1, 'r')
plot(x1,y1_med, 'g')
xlim([100 700])
ylim([100 700])
title('Quadratic Fit of Reduced Data (taking median vs. mean of an interval)')
legend('Reduced Data','Quadratic Fit (mean of interval)','Quadratic Fit (med of interval)')
caption = sprintf('y = %f * x^2 + %f * x + %f', p(1), p(2), p(3));
text(150, 640, caption, 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
caption2 = sprintf('Mean: (%d clicks, compare with %d original clicks, decrease of %f%%)', length(row1), length(row), 100*(1-length(row1)/length(row)));
text(150, 610, caption2, 'FontSize', 11, 'Color', 'black', 'FontWeight', 'bold');
caption3 = sprintf('y = %f * x^2 + %f * x + %f', p_med(1), p_med(2), p_med(3));
text(150, 570, caption3, 'FontSize', 11, 'Color', 'g', 'FontWeight', 'bold');
caption4 = sprintf('Median: (%d clicks, compare with %d original clicks, decrease of %f%%)', length(row1med), length(row), 100*(1-length(row1med)/length(row)));
text(150, 540, caption4, 'FontSize', 11, 'Color', 'black', 'FontWeight', 'bold');

% Reduced data - %keep points which are 1 standard deviation from MEAN
figure
plot(row1,col1,'b.')
hold on
plot(x1,y1, 'r')
plot(x1,y1+2*delta,'m--',x1,y1-2*delta,'m--')
xlim([100 700])
ylim([100 700])
title('Quadratic Fit of Reduced Data with 95% Prediction Interval')
legend('Reduced Data','Quadratic Fit','95% Prediction Interval')
caption = sprintf('y = %f * x^2 + %f * x + %f', p(1), p(2), p(3));
text(150, 650, caption, 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
caption2 = sprintf('(%d clicks, compare with %d original clicks, decrease of %f%%)', length(row1), length(row), 100*(1-length(row1)/length(row)));


%% Part 2: Fit all participant's responses to a polynomial, save y-values to run KS test in R

%create matrix of y_values
participants_y = zeros(r, 800); %choose any number for y's
x_participant = linspace(200,600,800);

for i = 1:r
    [row_part, col_part] = find(matrices(:,:,i));
    p_participant = polyfit(row_part,col_part,2);
    participants_y(i,:) = p_participant(1)*x_participant.^2 + p_participant(2)*x_participant + p_participant(3);
end

%each row is a participant, columns are y-values for KS test in R
xlswrite('y_values_of_participants_eye1',participants_y)

%% Part 3: Plot expert's clicks and fit a quadratic
eye = readtable('Upper eyelid excel data points/1_1.xlsx');
%change the table to a matrix
eye_coord_mat = round(table2array(eye));

%create a matrix
eye_mat = zeros(1000,1000);

for i = 1:length(eye_coord_mat)
    x = eye_coord_mat(i,1);
    y = eye_coord_mat(i,2);
    eye_mat(x,y)= 1;
end

I_eye = mat2gray(transpose(eye_mat));
imshow(I_eye)

%fit a curve to expert

[row, col] = find(eye_mat);
[p_expert, S] = polyfit(row,col,2);
x1_exp = linspace(200, 600, 800);
[y1_exp, delta_exp] = polyval(p_expert, x1_exp, S);

figure
plot(row,col,'b.')
hold on
plot(x1_exp,y1_exp, 'r')
plot(x1_exp,y1_exp+2*delta_exp,'m--',x1_exp,y1_exp-2*delta_exp,'m--')
xlim([100 700])
ylim([100 700])
title('Quadratic Fit of Expert Data with 95% Prediction Interval')
legend('Expert Data','Quadratic Fit','95% Prediction Interval')
caption = sprintf('y = %f * x^2 + %f * x + %f', p_expert(1), p_expert(2), p_expert(3));
text(150, 650, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');

hold off

% Plot expert and participant 1
%x values from 200 to 600, try 802 values
xrange = linspace(200, 600, 800);

%expected y from gold standard
y_exp = p_expert(1)*xrange.^2 + p_expert(2)*xrange + p_expert(3);
%observed y 
y_obs = p(1)*xrange.^2 + p(2)*xrange + p(3);

plot(xrange, y_exp, '-r');
hold on
plot(xrange, y_obs, '-b')
xlim([100 700])
ylim([100 700])
title('Quadratic Fit of Expert Data and Average Data');
legend( 'Expert fit', 'Average fit');
caption = sprintf('y = %f * x^2 + %f * x + %f', p_expert(1), p_expert(2), p_expert(3));
text(150, 660, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
caption = sprintf('y = %f * x^2 + %f * x + %f', p(1), p(2), p(3));
text(150, 620, caption, 'FontSize', 14, 'Color', 'b', 'FontWeight', 'bold');
hold off

%% Part 4: Generate a sample participant's response vs expert
[row, col] = find(matrices(:,:,33));
[p_individual, S] = polyfit(row,col,2);
x1 = linspace(200,600, 800);
[y1, delta] = polyval(p_individual, x1, S);

figure
plot(row,col,'b.')
hold on
plot(x1,y1, 'b')
plot(x1_exp, y_exp, '-r');
plot(x1_exp,y1_exp+2*delta_exp,'m--',x1_exp,y1_exp-2*delta_exp,'m--')
xlim([100 700])
ylim([100 700])
%rectangle('position',[200 100 10 600])
title('Quadratic Fit of Data with 95% Prediction Interval')
legend('Participant Response','Participant Fit', 'Expert Fit','95% Expert PI')
caption = sprintf('y = %f * x^2 + %f * x + %f', p_individual(1), p_individual(2), p_individual(3));
text(150, 650, caption, 'FontSize', 14, 'Color', 'b', 'FontWeight', 'bold');
caption = sprintf('y = %f * x^2 + %f * x + %f', p_expert(1), p_expert(2), p_expert(3));
text(150, 620, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');

hold off

%observed y 
y_obs_single = p_individual(1)*xrange.^2 + p_individual(2)*xrange + p_individual(3);

plot(xrange, y_exp, '-r');
hold on
plot(xrange, y_obs_single)
xlim([100 700])
ylim([100 700])
title('Quadratic Fit of Expert Data and Participant Data');
legend( 'Expert fit', 'Participant fit');
hold off

%% Part 5: Analysis of Number of Clicks

num_clicks = zeros(r, 1);
%calculate number of clicks
for i = 1:r
    temp = find(rootman_mat(i,:));
    num_clicks(i) = length(temp)/2;
end

%find mean number of clicks
quart = quantile(num_clicks,[0.25, 0.5, 0.75])

mean_clicks = round(mean(num_clicks));

%find the indices of participants which fall below and above the average
%number of clicks
below_mean_indx = find(num_clicks < mean_clicks);
above_mean_indx = find(num_clicks >= mean_clicks);

below_mean_indx = find(num_clicks < mean_clicks);
above_mean_indx = find(num_clicks >= mean_clicks);

%create two knew matrices for both groups
below_mean = zeros(length(below_mean_indx), c);
above_mean = zeros(length(above_mean_indx), c);

%fill in the matrices
for i = 1:length(below_mean_indx)
    below_mean(i, :) = rootman_mat(below_mean_indx(i), :);
end

for i = 1:length(above_mean_indx)
    above_mean(i, :) = rootman_mat(above_mean_indx(i), :);
end

%for all of the participants below the mean of clicks create a 3D matrix
mat_below = zeros(1000,1000,length(below_mean_indx));
for i = 1:length(below_mean_indx) %go through EACH participant
    for j = 1:(c-1)
        if mod(j , 2) == 1
           x = below_mean(i,j);
           j = j+1;
           y = below_mean(i,j);
           if x ~= 0 && y~= 0
              mat_below(x,y,i)= 1;
           end       
        end   
    end
end

mat_below_summed = (sum(mat_below,3))/length(below_mean_indx);

% for all of the participants avove the mean of clicks create a matrix
mat_above = zeros(1000,1000,length(above_mean_indx));
for i = 1:length(above_mean_indx)
    for j = 1:(c-1)
        if mod(j , 2) == 1
           x = above_mean(i,j);
           j = j+1;
           y = above_mean(i,j);
           if x ~= 0 && y~= 0
              mat_above(x,y,i)= 1;
           end       
        end   
    end
end

mat_above_summed = (sum(mat_above,3))/length(above_mean_indx);

%fit a curve to below average group

x1 = linspace(200,600,800);

[row_below, col_below] = find(mat_below_summed);
[p_below, S_below] = polyfit(row_below,col_below,2);
[y_below, delta] = polyval(p_below, x1, S_below);

figure
plot(row_below,col_below,'b.')
hold on
plot(x1,y_below, 'r')
plot(x1,y_below+2*delta,'m--',x1,y_below-2*delta,'m--')
xlim([0 800])
ylim([0 800])
%rectangle('position',[200 100 10 600])
title('Quadratic Fit of Averaged Participant Response (less than Average # of Clicks)')
legend('Data','Quadratic Fit','95% Prediction Interval')
caption = sprintf('y = %f * x^2 + %f * x + %f', p_below(1), p_below(2), p_below(3));
text(50, 750, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');

hold off

%fit a curve to REDUCED below average group
std_colmns = zeros(44,1);
j = 1;
for i = [181 191 201 211 221 231 241 251 261 271 281 291 301 311 321 331 341 351 361 371 381 391 401 411 421 431 441 451 461 471 481 491 501 511 521 531 541 551 561 571 581 591 601 611]
    [xx, yy] = find(mat_below_summed(i:(i+9),:));
    std_colmns(j) = std(yy);
    j=j+1;
end

std_colmns2 = std_colmns(~isnan(std_colmns));
mean_std = mean(std_colmns2);

%try mean
mat_below_summed2 = zeros(1000,1000);

mean_yy = 0;
minusSD_yy = 0;
plusSD_yy = 0;

for i = [181 191 201 211 221 231 241 251 261 271 281 291 301 311 321 331 341 351 361 371 381 391 401 411 421 431 441 451 461 471 481 491 501 511 521 531 541 551 561 571 581 591 601 611]
    [xx, yy, z] = find(mat_below_summed(i:(i+9),:));
    mean_yy = mean(yy);
    minusSD_yy = mean_yy - mean_std;
    plusSD_yy = mean_yy + mean_std;
    for j = 1:length(yy)
        if (yy(j) < plusSD_yy) && (yy(j) > minusSD_yy)
            mat_below_summed2((xx(j) + i - 1), yy(j)) = z(j);
        end
    end    
end

[row_below2, col_below2] = find(mat_below_summed2);
[p_below2, S_below2] = polyfit(row_below2,col_below2,2);
[y_below2, delta2] = polyval(p_below2, x1, S_below2);

figure
plot(row_below2,col_below2,'b.')
hold on
plot(x1,y_below2, 'r')
plot(x1,y_below2+2*delta2,'m--',x1,y_below2-2*delta2,'m--')
xlim([0 800])
ylim([0 800])
%rectangle('position',[200 100 10 600])
title('Quadratic Fit of REDUCED Averaged Participant Response (less than Average # of Clicks)')
legend('Reduced Data','Quadratic Fit','95% Prediction Interval')
caption = sprintf('y = %f * x^2 + %f * x + %f', p_below2(1), p_below2(2), p_below2(3));
text(50, 750, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
caption2 = sprintf('(%d clicks, compare with %d original clicks, decrease of %f%%)', length(row_below2), length(row_below), 100*(1-length(row_below2)/length(row_below)));
text(50, 700, caption2, 'FontSize', 11, 'Color', 'black', 'FontWeight', 'bold');

hold off

%fit a curve to above average group

[row_above, col_above] = find(mat_above_summed);
[p_above, S_above] = polyfit(row_above,col_above,2);
[y_above, delta] = polyval(p_above, x1, S_above);

figure
plot(row_above,col_above,'b.')
hold on
plot(x1,y_above, 'r')
plot(x1,y_above+2*delta,'m--',x1,y_above-2*delta,'m--')
xlim([0 800])
ylim([0 800])
%rectangle('position',[200 100 10 600])
title('Quadratic Fit of Averaged Participant Response (more than Average # of Clicks)')
legend('Data','Quadratic Fit','95% Prediction Interval')
caption = sprintf('y = %f * x^2 + %f * x + %f', p_above(1), p_above(2), p_above(3));
text(50, 750, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');

hold off

%fit a curve to REDUCED above average group
std_colmns = zeros(44,1);
j = 1;
for i = [181 191 201 211 221 231 241 251 261 271 281 291 301 311 321 331 341 351 361 371 381 391 401 411 421 431 441 451 461 471 481 491 501 511 521 531 541 551 561 571 581 591 601 611]
    [xx, yy] = find(mat_above_summed(i:(i+9),:));
    std_colmns(j) = std(yy);
    j=j+1;
end

std_colmns2 = std_colmns(~isnan(std_colmns));
mean_std = mean(std_colmns2);

%try mean
mat_above_summed2 = zeros(1000,1000);

mean_yy = 0;
minusSD_yy = 0;
plusSD_yy = 0;

for i = [181 191 201 211 221 231 241 251 261 271 281 291 301 311 321 331 341 351 361 371 381 391 401 411 421 431 441 451 461 471 481 491 501 511 521 531 541 551 561 571 581 591 601 611]
    [xx, yy, z] = find(mat_above_summed(i:(i+9),:));
    mean_yy = mean(yy);
    minusSD_yy = mean_yy - mean_std;
    plusSD_yy = mean_yy + mean_std;
    for j = 1:length(yy)
        if (yy(j) < plusSD_yy) && (yy(j) > minusSD_yy)
            mat_above_summed2((xx(j) + i - 1), yy(j)) = z(j);
        end
    end    
end

[row_above2, col_above2] = find(mat_above_summed2);
[p_above2, S_above2] = polyfit(row_above2,col_above2,2);
[y_above2, delta2] = polyval(p_above2, x1, S_above2);

figure
plot(row_above2,col_above2,'b.')
hold on
plot(x1,y_above2, 'r')
plot(x1,y_above2+2*delta2,'m--',x1,y_above2-2*delta2,'m--')
xlim([0 800])
ylim([0 800])
%rectangle('position',[200 100 10 600])
title('Quadratic Fit of REDUCED Averaged Participant Response (more than Average # of Clicks)')
legend('Reduced Data','Quadratic Fit','95% Prediction Interval')
caption = sprintf('y = %f * x^2 + %f * x + %f', p_above2(1), p_above2(2), p_above2(3));
text(50, 750, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
caption2 = sprintf('(%d clicks, compare with %d original clicks, decrease of %f%%)', length(row_above2), length(row_above), 100*(1-length(row_above2)/length(row_above)));
text(50, 700, caption2, 'FontSize', 11, 'Color', 'black', 'FontWeight', 'bold');

hold off

%Plot all on same graph
y_exp2 = p_expert(1)*x1.^2 + p_expert(2)*x1 + p_expert(3);

figure
plot(x1,y_below2,'r')
hold on
plot(x1,y_above2, 'b')
plot(x1,y_exp2,'g')
xlim([100 700])
ylim([100 700])
title('Quadratic Fit of REDUCED Avg Participant Response')
%legend('Below Avg Reduced Fit','Above Avg Reduced Fit','Expert Fit')
legend('(Less than 10 clicks) Avg Reduced Fit','(15 or more clicks) Avg Reduced Fit','Expert Fit')
caption = sprintf('y = %f * x^2 + %f * x + %f', p_below2(1), p_below2(2), p_below2(3));
text(150, 650, caption, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');

caption2 = sprintf('y = %f * x^2 + %f * x + %f', p_above2(1), p_above2(2), p_above2(3));
text(150, 625, caption2, 'FontSize', 14, 'Color', 'b', 'FontWeight', 'bold');

caption3 = sprintf('y = %f * x^2 + %f * x + %f', p_expert(1), p_expert(2), p_expert(3));
text(150, 600, caption3, 'FontSize', 14, 'Color', 'g', 'FontWeight', 'bold');



