%% This code runs through the whole processing pipeline and facilitates a user input

% Created by: Dylan Beam on 3/23/2020
% Last Edited by: Dylan Beam on 3/23/2020

%% Load Metadata and define variables

clc,clear

%load metadata
load('expmt_list_1_6_20.mat')

%define the user path for the computer you're using
source_path = 'R:\users\dyb11\source';

%Define constants across animals
c.fs=30e3; %sampling frequency
c.time = 60; %the duration of a single trial (milliseconds)
c.bin = [-2 400]; %window for STA calculation where stim event is 0
c.t_min = c.bin(1)*c.fs/1000; %convert ms to s, then multiply by sample rate
c.t_max = c.bin(2)*c.fs/1000;
c.window_size = c.t_max - c.t_min; %the number of sampled points in the appropriate window size
c.buffer = 1.5; %time (ms) around the stimulation event to be removed

%% Segment data for processing

data_table = segment_raw_data(expmt_list, 3)

%% Collect baseline values

%loop through all animals
for expmt = 2:7
    
    %define animal specific parameters
    exp_data = expmt_list{expmt}; %load animal specific data
    stim_port = exp_data.stim_port; %identify the stim port for this animal
    cohort = sprintf('F%s', exp_data.cohort); %identification number of the ferret
    rec_path = [source_path, cohort, '\Trellis']; %access path to data
    num_trials = exp_data.trial_list(6,2); %number of trials to analyze
    FailedTrials = zeros(32,num_trials); %matrix of trials that have failed
    break_point = exp_data.trial_list(3,2); %the trial number where stimulation switches from 1:2 to 3:4
    switch expmt %identify trials to be excluded for the given animal
        case 2
            excluded_trials = [1, 29:33, 42:48];
        case 3
            excluded_trials = [31:37, 45:48];
        case 4
            excluded_trials = [1, 6, 24];
        case 5
            excluded_trials = [1, 16:18];
        case 6
            excluded_trials = [26:28];
        case 7
            excluded_trials = [1, 23];
    end
    
    %Collect baseline noise data for the designated animal and visualize
    %results
    [baseline] = visualizeNoise(source_path, expmt_list, expmt, excluded_trials);
        
    %finish collecting baseline data from failed trials
    finishing = 'y';
    while finishing == 'y'
        finishing = input('Were there failed trials that need to be reprocessed? \n');
        
        %[baseline] = finishBaseline(source_path, expmt_list, expmt, excluded_trials);
        
    end   
    
end

%% Detect Responses

%% Create Selectivity data

%% Create Selective Stimulation Boundaries

