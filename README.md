# EH-deploy

Reliability-Driven Deployment in Energy-Harvesting Sensor Networks.

## Getting Started

Test environment: MATLAB R2019b/R2020a.

**Note:** Need [[CPLEX]](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.studio.help/pdf/gscplexmatlab.pdf) for Matlab for the optimization solver.

The `tutorial.m` will walk you through all of our algorithms on the small dataset.

## File Structure

```
.
├── LICENSE
├── README.md     // this file
├── exp           // scripts to run experiments
├── solardata     // real-world dataset from the National Solar Radiation Database
├── libs          // general library
├── lldistkm      // distance calculation library from FileExchange
├── solver        // the CPLEX solver
├── RDTSH.m       // our purposed heuristic, Reliability-Driven Two-StageHeuristic
├── SRIGH.m       // Sensing- and Routing- Integrated Greedy Heuristic
├── TSH.m         // Two-Stage Heuristic
└── tutorial.m
```

## Available Algorithms

* Reliability-Driven Two-StageHeuristic (RDTSH).
* Sensing- and Routing- Integrated Greedy Heuristic (SRIGH).
* Two-Stage Heuristic (TSH).

For all of the three algorithms, we devise the power model for our sensor deployment problem in `./solver/power_eno.m`.

## License

MIT