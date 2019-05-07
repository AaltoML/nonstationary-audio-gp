function [y] = GammaToneApply(x,forward,feedback)
% [y] = GammaToneApply(x,forward,feedback)

% Stanley 1993 Edited June 2007 by Nick C.
% I DID NOT DERIVE ALL THE CLEVER STUFF!! THAT WONDERFUL STUFF CAN BE
% FOUND IN Apple TR #35, "An Efficient Implementation of the Patterson-
% Holdsworth Cochlear Filter Bank."

% This function filters the waveform x with the array of filters
% specified by the forward and feedback parameters. Each row
% of the forward and feedback parameters are the parameters
% to the Matlab builtin function "filter".

[rows, cols]=size(feedback);
y=zeros(rows,length(x));
for i=1:rows
    y(i,:)=filter(forward(i,:),feedback(i,:),x);
end

