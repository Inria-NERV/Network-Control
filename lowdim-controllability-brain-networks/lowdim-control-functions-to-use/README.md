# Ready-to-use MATLAB functions in order to characterize the control centrality of driver nodes in an undirected network.

The function "low_dimensional_control_centrality" answers the problem of selecting the best single drivers when there are no states specified but rather in a general sense (indeed we consider the worst-case direction).

The function "controlMetric_forStateTransfer" answers the same question of selecting drivers but for a specific case where the user specifies initial and final states x0 and xf.

--> Both functions return two metrics; the standard one, and its low-dimensional version (The method explained in the publication).

