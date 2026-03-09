% Data from Question 6
days = [0 6 10 13 17 20 28];
young = [7 18 45 40 32 30.5 30];
mature = [7 16 19 15 12 10.5 10];

%% (a)

xstar = 3;
young_est = cubic_spline(days,young,xstar);
mature_est = cubic_spline(days,mature,xstar);

% built-in spline (not-a-knot)
young_nak = spline(days,young,xstar);
mature_nak = spline(days,mature,xstar);

disp('Weight at 3 days (young leaves):')
disp(young_est)
disp('Weight at 3 days (young leaves, built-in spline):')
disp(young_nak)
disp('Weight at 3 days (mature leaves):')
disp(mature_est)
disp('Weight at 3 days (mature leaves, built-in spline):')
disp(mature_nak)

%% (b)

% create evaluation points
xvals = linspace(0,28,200);

young_vals = cubic_spline(days,young,xvals);
mature_vals = cubic_spline(days,mature,xvals);

[max_young, idx_young] = max(young_vals);
[max_mature, idx_mature] = max(mature_vals);

day_young = xvals(idx_young);
day_mature = xvals(idx_mature);

disp('Maximum weight (young leaves):')
disp(max_young)
disp('Occurs at day:')
disp(day_young)
disp('Maximum weight (mature leaves):')
disp(max_mature)
disp('Occurs at day:')
disp(day_mature)