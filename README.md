# research-pse
Partial Speculative Execution, abbrev. PSE

This project is to guide a implmentation of PSE on top of Hadoop-0.21.0. The motivation and experimental results will be linked to a research paper. 

PSE is designed to prevent speculative tasks from re-reading, re-copying, and re-computing the processed data, while improving the efficiency of speculative execution and decreasing job completion time.

1. Lightweight Check-pointing and Recovery Algorithm
2. State Transition //see the paper
3. The PSE Scheduler //see the paper
