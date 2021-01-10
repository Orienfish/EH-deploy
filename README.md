# EH-deploy

[CNSM 2020] Reliability-Driven Deployment in Energy-Harvesting Sensor Networks.

## Getting Started

Test environment: MATLAB R2019b/R2020a.

**Note:** Need [[CPLEX]](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.studio.help/pdf/gscplexmatlab.pdf) for MATLAB as the optimization solver.

The `tutorial.m` will walk you through all of our algorithms on the dataset.

## File Structure

```
.
├── LICENSE
├── README.md     // this file
├── exp           // scripts to run experiments
├── solardata     // real-world dataset from the National Solar Radiation Database (NSRDB)
├── libs          // general library for reliability models
├── lldistkm      // distance calculation library from FileExchange
├── solver        // the matlab script to call the CPLEX solver
├── alg			      // our purposed heuristic, Reliability-Driven Two-StageHeuristic
└── tutorial.m    // example usage
```

## Available Algorithms

Currently we have the following heuristics under the `./alg/` directory. 

* Reliability-Driven Two-Stage Heuristic (RDTSH) [our proposed heuristic]. 
* Two-Stage Heuristic (TSH) [[Zhu et al. 2018]](https://ieeexplore.ieee.org/abstract/document/8345168).
* Sensing- and Routing- Integrated Greedy Heuristic (SRIGH) [[Zhu et al. 2018]](https://ieeexplore.ieee.org/abstract/document/8345168).

## License

MIT

