%EEGLAB graphic interface functions (popfunc folder):
%  <a href="matlab:helpwin eeg_addnewevents">eeg_addnewevents</a>     - Eeg_addnewevents() Add new events to EEG structure. Both EEG.event and...
%  <a href="matlab:helpwin eeg_amplitudearea">eeg_amplitudearea</a>    - Resamples an ERP average using spline interpolation...
%  <a href="matlab:helpwin eeg_chaninds">eeg_chaninds</a>         - Look up channel indices in a EEG structure...
%  <a href="matlab:helpwin eeg_chantype">eeg_chantype</a>         - Returns the channel indices of the desired channel type(s).
%  <a href="matlab:helpwin eeg_context">eeg_context</a>          - Returns (in output 'delays') a matrix giving, for each event of specified...
%  <a href="matlab:helpwin eeg_decodechan">eeg_decodechan</a>       - Given an input EEG dataset structure, output a new EEG data structure...
%  <a href="matlab:helpwin eeg_dipselect">eeg_dipselect</a>        - Select componet dipoles from an EEG dataset with...
%  <a href="matlab:helpwin eeg_eegrej">eeg_eegrej</a>           - Reject porition of continuous data in an EEGLAB...
%  <a href="matlab:helpwin eeg_emptyset">eeg_emptyset</a>         - Initialize an EEG dataset structure with default values.
%  <a href="matlab:helpwin eeg_epoch2continuous">eeg_epoch2continuous</a> - Convert epoched dataset to continuous dataset...
%  <a href="matlab:helpwin eeg_epochformat">eeg_epochformat</a>      - Convert the epoch information of a dataset from struct...
%  <a href="matlab:helpwin eeg_eventformat">eeg_eventformat</a>      - Convert the event information of a dataset from struct...
%  <a href="matlab:helpwin eeg_eventhist">eeg_eventhist</a>        - Return or plot histogram of event or urevent field values.
%  <a href="matlab:helpwin eeg_eventtable">eeg_eventtable</a>       - Returns all events contained in the EEG structure (and...
%  <a href="matlab:helpwin eeg_eventtypes">eeg_eventtypes</a>       - Return a list of event or urevent types in a dataset and...
%  <a href="matlab:helpwin eeg_getepochevent">eeg_getepochevent</a>    - Return dataset event field values for all events...
%  <a href="matlab:helpwin eeg_getica">eeg_getica</a>           - Get ICA component activation. Recompute if necessary.
%  <a href="matlab:helpwin eeg_insertbound">eeg_insertbound</a>      - Insert boundary event in an EEG event structure.
%  <a href="matlab:helpwin eeg_interp">eeg_interp</a>           - Interpolate data channels...
%  <a href="matlab:helpwin eeg_laplac">eeg_laplac</a>           - Gives the laplacian for the data contained in EEG...
%  <a href="matlab:helpwin eeg_lat2point">eeg_lat2point</a>        - Convert latencies in time units relative to the...
%  <a href="matlab:helpwin eeg_latencyur">eeg_latencyur</a>        - Transform latency of sample point in the continuous...
%  <a href="matlab:helpwin eeg_matchchans">eeg_matchchans</a>       - Find closest channels in a larger EEGLAB chanlocs structure...
%  <a href="matlab:helpwin eeg_mergechan">eeg_mergechan</a>        - Merge channel structure while preserving channel...
%  <a href="matlab:helpwin eeg_mergelocs">eeg_mergelocs</a>        - Merge channel structure while preserving channel...
%  <a href="matlab:helpwin eeg_mergelocs_diffstruct">eeg_mergelocs_diffstruct</a> - Merge channel structure while preserving channel...
%  <a href="matlab:helpwin eeg_multieegplot">eeg_multieegplot</a>     - Produce an eegplot() of a the average of an epoched dataset...
%  <a href="matlab:helpwin eeg_oldica">eeg_oldica</a>           - Report, return or add to oldicaweights and oldicasphere...
%  <a href="matlab:helpwin eeg_point2lat">eeg_point2lat</a>        - Convert latency in data points to latency in ms relative...
%  <a href="matlab:helpwin eeg_pv">eeg_pv</a>               - Compute EEG.data 'percent variance ' (pv) of the whole EEG data versus the projections...
%  <a href="matlab:helpwin eeg_pvaf">eeg_pvaf</a>             - Compute EEG.data 'percent variance accounted for' (pvaf) by specified components.
%  <a href="matlab:helpwin eeg_rejmacro">eeg_rejmacro</a>         - Internal EEGLAB macro for all pop_ functions that...
%  <a href="matlab:helpwin eeg_rejsuperpose">eeg_rejsuperpose</a>     - Superpose rejections of a EEG dataset.
%  <a href="matlab:helpwin eeg_timeinterp">eeg_timeinterp</a>       - Perform spline interpolation of a portion...
%  <a href="matlab:helpwin eeg_topoplot">eeg_topoplot</a>         - Plot scalp map...
%  <a href="matlab:helpwin eeg_urlatency">eeg_urlatency</a>        - Find the original (ur) latency of a time point in...
%  <a href="matlab:helpwin getchanlist">getchanlist</a>          - Obtain indices of specified channel types.
%  <a href="matlab:helpwin importevent">importevent</a>          - Import experimental events from data file or Matlab...
%  <a href="matlab:helpwin pop_autorej">pop_autorej</a>          - Perform automatic artifact epoch detection and rejection...
%  <a href="matlab:helpwin pop_averef">pop_averef</a>           - Convert an EEG dataset to average reference.
%  <a href="matlab:helpwin pop_biosig">pop_biosig</a>           - Import data files into EEGLAB using BIOSIG toolbox...
%  <a href="matlab:helpwin pop_biosig16">pop_biosig16</a>         - Import data files into EEGLAB using BIOSIG toolbox...
%  <a href="matlab:helpwin pop_biosig16ying">pop_biosig16ying</a>     - Import data files into EEGLAB using BIOSIG toolbox...
%  <a href="matlab:helpwin pop_chancenter">pop_chancenter</a>       - Recenter cartesian X,Y,Z channel coordinates...
%  <a href="matlab:helpwin pop_chancoresp">pop_chancoresp</a>       - Define correspondances between two channel locations structures...
%  <a href="matlab:helpwin pop_chanedit">pop_chanedit</a>         - Edit the channel locations structure of an EEGLAB dataset,...
%  <a href="matlab:helpwin pop_chanevent">pop_chanevent</a>        - Import event latencies from the rising and/or falling 'edge'...
%  <a href="matlab:helpwin pop_chansel">pop_chansel</a>          - Pop up a graphic interface to select channels...
%  <a href="matlab:helpwin pop_comments">pop_comments</a>         - Edit comments...
%  <a href="matlab:helpwin pop_compareerps">pop_compareerps</a>      - Compare the (ERP) averages of two datasets.
%  <a href="matlab:helpwin pop_comperp">pop_comperp</a>          - Compute the grand average ERP waveforms of multiple datasets...
%  <a href="matlab:helpwin pop_copyset">pop_copyset</a>          - Copy the current EEG dataset into another dataset.
%  <a href="matlab:helpwin pop_crossf">pop_crossf</a>           - Return estimates and plots of event-related spectral coherence...
%  <a href="matlab:helpwin pop_editeventfield">pop_editeventfield</a>   - Add/remove/rename/modify a field in the event structure...
%  <a href="matlab:helpwin pop_editeventvals">pop_editeventvals</a>    - Edit events contained in an EEG dataset structure.
%  <a href="matlab:helpwin pop_editset">pop_editset</a>          - Edit EEG dataset structure fields.
%  <a href="matlab:helpwin pop_eegfilt">pop_eegfilt</a>          - Interactively filter EEG dataset data using eegfilt()...
%  <a href="matlab:helpwin pop_eegplot">pop_eegplot</a>          - Visually inspect EEG data using a scrolling display.
%  <a href="matlab:helpwin pop_eegthresh">pop_eegthresh</a>        - Reject artifacts by detecting outlier values.  This has...
%  <a href="matlab:helpwin pop_envtopo">pop_envtopo</a>          - Plot envelope of an averaged EEG epoch, plus scalp maps...
%  <a href="matlab:helpwin pop_epoch">pop_epoch</a>            - Convert a continuous EEG dataset to epoched data by extracting...
%  <a href="matlab:helpwin pop_erpimage">pop_erpimage</a>         - Draw an ERP-image plot of a given EEG channel or independent...
%  <a href="matlab:helpwin pop_eventstat">pop_eventstat</a>        - Computes and plots statistical characteristics of an EEG event,...
%  <a href="matlab:helpwin pop_expevents">pop_expevents</a>        - Export events to CSV file...
%  <a href="matlab:helpwin pop_expica">pop_expica</a>           - Export ICA weights or inverse matrix...
%  <a href="matlab:helpwin pop_export">pop_export</a>           - Export EEG dataset...
%  <a href="matlab:helpwin pop_fileio">pop_fileio</a>           - Import data files into EEGLAB using FileIO...
%  <a href="matlab:helpwin pop_fileiodir">pop_fileiodir</a>        - Import directory into EEGLAB using FileIO...
%  <a href="matlab:helpwin pop_headplot">pop_headplot</a>         - Plot one or more spherically-splined EEG field maps...
%  <a href="matlab:helpwin pop_icathresh">pop_icathresh</a>        - Main menu for choosing threshold for component...
%  <a href="matlab:helpwin pop_importdata">pop_importdata</a>       - Import data from a Matlab variable or disk file by calling...
%  <a href="matlab:helpwin pop_importegimat">pop_importegimat</a>     - Import EGI Matlab segmented file...
%  <a href="matlab:helpwin pop_importepoch">pop_importepoch</a>      - Export epoch and/or epoch event information to the event...
%  <a href="matlab:helpwin pop_importev2">pop_importev2</a>        - Merge a neuroscan EV2 file with input dataset...
%  <a href="matlab:helpwin pop_importevent">pop_importevent</a>      - Import events into an EEG dataset. If the EEG dataset...
%  <a href="matlab:helpwin pop_importpres">pop_importpres</a>       - Append Presentation event file information into an EEGLAB dataset...
%  <a href="matlab:helpwin pop_interp">pop_interp</a>           - Interpolate data channels...
%  <a href="matlab:helpwin pop_jointprob">pop_jointprob</a>        - Reject artifacts in an EEG dataset using joint...
%  <a href="matlab:helpwin pop_loadbci">pop_loadbci</a>          - Import a BCI2000 ascii file into EEGLAB...
%  <a href="matlab:helpwin pop_loadcnt">pop_loadcnt</a>          - Load a neuroscan CNT file (pop out window if no arguments).
%  <a href="matlab:helpwin pop_loaddat">pop_loaddat</a>          - Merge a neuroscan DAT file with input dataset...
%  <a href="matlab:helpwin pop_loadeeg">pop_loadeeg</a>          - Load a Neuroscan .EEG file (via a pop-up window if no...
%  <a href="matlab:helpwin pop_loadset">pop_loadset</a>          - Load an EEG dataset. If no arguments, pop up an input window.
%  <a href="matlab:helpwin pop_mergeset">pop_mergeset</a>         - Merge two or more datasets. If only one argument is given,...
%  <a href="matlab:helpwin pop_newcrossf">pop_newcrossf</a>        - Return estimates and plots of event-related spectral coherence...
%  <a href="matlab:helpwin pop_newset">pop_newset</a>           - Edit/save EEG dataset structure information.
%  <a href="matlab:helpwin pop_newtimef">pop_newtimef</a>         - Returns estimates and plots of event-related (log) spectral...
%  <a href="matlab:helpwin pop_plotdata">pop_plotdata</a>         - Plot average of EEG channels or independent components in...
%  <a href="matlab:helpwin pop_plottopo">pop_plottopo</a>         - Plot one or more concatenated multichannel data epochs...
%  <a href="matlab:helpwin pop_prop">pop_prop</a>             - Plot the properties of a channel or of an independent...
%  <a href="matlab:helpwin pop_readegi">pop_readegi</a>          - Load a EGI EEG file (pop out window if no arguments).
%  <a href="matlab:helpwin pop_readlocs">pop_readlocs</a>         - Load a EGI-format EEG file (pop up an interactive window if no arguments).
%  <a href="matlab:helpwin pop_readsegegi">pop_readsegegi</a>       - Load a segmented EGI EEG file. Pop up query...
%  <a href="matlab:helpwin pop_rejchan">pop_rejchan</a>          - Reject artifacts channels in an EEG dataset using joint...
%  <a href="matlab:helpwin pop_rejchanspec">pop_rejchanspec</a>      - Reject artifacts channels in an EEG dataset using...
%  <a href="matlab:helpwin pop_rejcont">pop_rejcont</a>          - Reject continuous portions of data based on spectrum...
%  <a href="matlab:helpwin pop_rejepoch">pop_rejepoch</a>         - Reject pre-labeled trials in a EEG dataset.
%  <a href="matlab:helpwin pop_rejkurt">pop_rejkurt</a>          - Rejection of artifact in a dataset using kurtosis...
%  <a href="matlab:helpwin pop_rejspec">pop_rejspec</a>          - Rejection of artifact in a dataset using...
%  <a href="matlab:helpwin pop_rejtrend">pop_rejtrend</a>         - Measure linear trends in EEG data; reject data epochs...
%  <a href="matlab:helpwin pop_reref">pop_reref</a>            - Convert an EEG dataset to average reference or to a...
%  <a href="matlab:helpwin pop_resample">pop_resample</a>         - Resample dataset (pop up window).
%  <a href="matlab:helpwin pop_rmbase">pop_rmbase</a>           - Remove channel baseline means from an epoched or...
%  <a href="matlab:helpwin pop_rmdat">pop_rmdat</a>            - Remove continuous data around specific events...
%  <a href="matlab:helpwin pop_runica">pop_runica</a>           - Run an ICA decomposition of an EEG dataset using runica(),...
%  <a href="matlab:helpwin pop_runscript">pop_runscript</a>        - Run Matlab script...
%  <a href="matlab:helpwin pop_saveh">pop_saveh</a>            - Save the EEGLAB session command history stored in ALLCOM...
%  <a href="matlab:helpwin pop_saveset">pop_saveset</a>          - Save one or more EEG dataset structures...
%  <a href="matlab:helpwin pop_select">pop_select</a>           - Given an input EEG dataset structure, output a new EEG data structure...
%  <a href="matlab:helpwin pop_selectcomps">pop_selectcomps</a>      - Display components with button to vizualize their...
%  <a href="matlab:helpwin pop_selectevent">pop_selectevent</a>      - Find events in an EEG dataset. If the dataset...
%  <a href="matlab:helpwin pop_signalstat">pop_signalstat</a>       - Computes and plots statistical characteristics of a signal,...
%  <a href="matlab:helpwin pop_snapread">pop_snapread</a>         - Load an EEG SnapMaster file (pop out window if no arguments).
%  <a href="matlab:helpwin pop_spectopo">pop_spectopo</a>         - Plot spectra of specified data channels or components.
%  <a href="matlab:helpwin pop_subcomp">pop_subcomp</a>          - Remove specified components from an EEG dataset.
%  <a href="matlab:helpwin pop_timef">pop_timef</a>            - Returns estimates and plots of event-related (log) spectral...
%  <a href="matlab:helpwin pop_timtopo">pop_timtopo</a>          - Call the timtopo() function for epoched EEG datasets.
%  <a href="matlab:helpwin pop_topoplot">pop_topoplot</a>         - Plot scalp map(s) in a figure window. If number of input...
%  <a href="matlab:helpwin pop_writeeeg">pop_writeeeg</a>         - Write EEGLAB dataset to disk in EDF/GDF or BDF format...
%  <a href="matlab:helpwin pop_writelocs">pop_writelocs</a>        - Load a EGI EEG file (pop out window if no arguments).
