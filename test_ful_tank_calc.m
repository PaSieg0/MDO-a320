%% Extended Tests for Wing Tank Volume Calculation
clc; clear; close all;

% --- 1. SETUP PARAMETERS ---
airfoil_file = 'b737a.dat';

% Tank Limits (% of chord)
limit_start = 0.15;
limit_end   = 0.85;

% Base Wing Geometry
c_r  = 6.0;
c_k  = 3.8;
c_t  = 1.5;
b_k  = 5.5;
b    = 34.1;

fprintf('=== WING TANK VOLUME TEST SUITE ===\n\n');

%% TEST 1: Base Case
fprintf('TEST 1: Base Case\n');
[Vol_Total, Vol_Inner, Vol_Outer] = calculate_wing_tank_volume( ...
    c_r, c_k, c_t, b_k, b, airfoil_file, limit_start, limit_end);
fprintf('Total Volume: %.4f m^3\n\n', Vol_Total);

%% TEST 2: Proportionality - Double all chords (Volume should scale ~x8)
fprintf('TEST 2: Proportionality - Double all chords\n');
scale = 2.0;
[Vol_Scaled, ~, ~] = calculate_wing_tank_volume( ...
    c_r*scale, c_k*scale, c_t*scale, b_k, b, airfoil_file, limit_start, limit_end);
ratio = Vol_Scaled / Vol_Total;
fprintf('Expected ratio: ~%.1f, Actual ratio: %.2f\n', scale^2, ratio);
assert(abs(ratio - scale^2) < 0.5, 'Chord scaling test FAILED');
fprintf('TEST 2: PASSED\n\n');

%% TEST 3: Proportionality - Double wingspan (Volume should scale ~x2)
fprintf('TEST 3: Proportionality - Double wingspan\n');
[Vol_SpanScaled, ~, ~] = calculate_wing_tank_volume( ...
    c_r, c_k, c_t, b_k*2, b*2, airfoil_file, limit_start, limit_end);
span_ratio = Vol_SpanScaled / Vol_Total;
fprintf('Expected ratio: ~2.0, Actual ratio: %.2f\n', span_ratio);
assert(abs(span_ratio - 2.0) < 0.3, 'Span scaling test FAILED');
fprintf('TEST 3: PASSED\n\n');

%% TEST 4: Tank limits - Wider limits should give more volume
fprintf('TEST 4: Tank Limits - Wider limits give more volume\n');
[Vol_Narrow, ~, ~] = calculate_wing_tank_volume( ...
    c_r, c_k, c_t, b_k, b, airfoil_file, 0.25, 0.75);
[Vol_Wide, ~, ~] = calculate_wing_tank_volume( ...
    c_r, c_k, c_t, b_k, b, airfoil_file, 0.10, 0.90);
fprintf('Narrow (25%%-75%%): %.4f m^3\n', Vol_Narrow);
fprintf('Wide (10%%-90%%):   %.4f m^3\n', Vol_Wide);
assert(Vol_Wide > Vol_Narrow, 'Tank limits test FAILED');
assert(Vol_Narrow < Vol_Total, 'Narrow limits should be less than base');
fprintf('TEST 4: PASSED\n\n');

%% TEST 5: Inner + Outer = Total
fprintf('TEST 5: Volume Consistency (Inner + Outer = Total)\n');
vol_sum = Vol_Inner + Vol_Outer;
diff_pct = abs(vol_sum - Vol_Total) / Vol_Total * 100;
fprintf('Inner: %.4f, Outer: %.4f, Sum: %.4f, Total: %.4f\n', ...
    Vol_Inner, Vol_Outer, vol_sum, Vol_Total);
fprintf('Difference: %.4f%%\n', diff_pct);
assert(diff_pct < 0.1, 'Volume consistency test FAILED');
fprintf('TEST 5: PASSED\n\n');

%% TEST 6: Non-negative volumes
fprintf('TEST 6: Non-negative Volumes\n');
assert(Vol_Total >= 0, 'Total volume should be non-negative');
assert(Vol_Inner >= 0, 'Inner volume should be non-negative');
assert(Vol_Outer >= 0, 'Outer volume should be non-negative');
fprintf('TEST 6: PASSED\n\n');

%% TEST 7: Symmetric scaling (half chord, half span -> 1/4 volume)
fprintf('TEST 7: Symmetric Scaling (0.5x chord, 0.5x span)\n');
[Vol_Half, ~, ~] = calculate_wing_tank_volume( ...
    c_r*0.5, c_k*0.5, c_t*0.5, b_k*0.5, b*0.5, airfoil_file, limit_start, limit_end);
expected_ratio = 0.5^3; % chord^2 * span
actual_ratio = Vol_Half / Vol_Total;
fprintf('Expected ratio: %.4f, Actual ratio: %.4f\n', expected_ratio, actual_ratio);
assert(abs(actual_ratio - expected_ratio) < 0.05, 'Symmetric scaling test FAILED');
fprintf('TEST 7: PASSED\n\n');

%% SUMMARY
fprintf('==============================\n');
fprintf('   ALL TESTS PASSED!\n');
fprintf('==============================\n');