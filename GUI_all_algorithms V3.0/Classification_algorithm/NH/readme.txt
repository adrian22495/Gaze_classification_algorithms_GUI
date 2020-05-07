README-file of the event detection algorithm accompanying the article
Nyström, M., & Holmqvist, K. (2010). An adaptive algorithm for fixation, saccade, and glissade detection in eyetracking data. Behavior research methods, 42(1), 188-204.

The algorithm is implemented in Matlab and tested on version R2017a.
Usage: Run the file: mainRun.m. Update the parameters in mainRun.m to fit your experiment.

Update 2018-06-11: Small corrections in the code have been made based on the comments in Friedman et al. (2018).

Friedman, L., Rigas, I., Abdulin, E., & Komogortsev, O. V. (2018). A novel evaluation of two related and two independent algorithms for eye movement classification during reading. Behavior research methods, 1-24.

If you have any questions or comments on the paper or the code, send me an email at
marcus.nystrom@humlab.lu.se. 

Observe that the algorithm is suitable ONLY for data collected from viewers
keeping their heads relatively still while watching static stimuli. It is not 
designed to handle data containing smooth pursuit movements. It was tested on 1250 Hz data from an SMI Hi-speed eye tracker (example data provided).