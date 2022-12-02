In Matlab R2022a, gifti.m, specifying a class, apparently can no longer be treated as a function, to load a given input file; MATLAB error: Execution of script gifti as a function is not supported. I use a workaround in hb_prune_graph.m, using the code function that is called within gifti.m to read the file. 

The two read_* were taken from spm12's gifti class private folder. Currently I just use: read_gifti_file.m, but the other one can also become handy if I come across a need to just directly read FreeSurfer's surface files.  

- xml_parser_hb.m same as xml_parser.m in spm12; just renamed.

- read_freesurfer_file_hb.m same as read_freesurfer_file.m in spm12; just renamed.

- There are some minimal changes in read_gifti_file_hb.m compared to read_gifti_file.m. 

[07.05.2022]

  

