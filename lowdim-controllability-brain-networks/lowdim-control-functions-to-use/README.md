The function "low_dimensional_control_centrality" answers the problem of selecting the best single drivers when there no states specified but rather in a general sense (indeed we cansider the worst-case direction).

The function "controlMetric_forStateTransfer" however, answers the same question of selecting drivers but for a specific case where the user specifies initial and final states x0 and xf.

--> both functions return two metrics; the standard one, and its low-dimensional version (The method explained in the publication).

