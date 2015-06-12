# research-pse
Partial Speculative Execution, abbrev. PSE

This project is to guide a implmentation of PSE on top of Hadoop-0.21.0. The motivation and experimental results will be linked to a research paper later. 

In the standard MapReduce framework, speculative tasks start from the ground up to do the same work as their original tasks. This from-the-ground-up execution which we name it as full speculative execution (FSE) has several disadvantages. First, a speculative map task has to re-read the input data, which increases the I/O cost since the input data always comes from Hadoop Distributed File System (HDFS). Second, a speculative reduce task needs to re-copy the intermediate data from almost all mappers, occupying extra network bandwidth. Third, all speculative tasks have to re-compute the partial processed data from disk and remote mappers. These re-reading, re-copying, and re-computing costs of speculative tasks make speculative execution in-efficient and do few contributions to job's completion time as expected.

PSE is designed to prevent speculative tasks from re-reading, re-copying, and re-computing the processed data, while improving the efficiency of speculative execution and decreasing job's completion time.

It contains three parts:
1. Lightweight Check-pointing and Recovery Algorithm
2. State Transition
3. The PSE Scheduler

Resources description:
1. conf folder: xml configuration file with new parameters.
2. org/hadoop/: critical modified java classes, based on hadoop-0.21.0, also contains a new folder 'pse', newly created classes.
3. pseSimulation folder: a simulation tool, to estimate the schedulers in terms of differect cluster resources.
