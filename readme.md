# Geodesic Graph
The data structure implemented in `source/geodesic_graph.h`, together with the visist

# Results
![](data/results.png)
Comparison between this geodesic solver, and the solver from [Crane et al. 2013](https://arxiv.org/abs/1204.6216), the fastest solver for approximate geodesiscs. For the latter, implementation provided by the author was used, using Cholmod as backend. Columns build report pre-processing times to pre-factor the system and to build the graph, respectively. Build times for this solver also include the time to compute the triangle adjacency needed to build the graph. Columns solve report average time for computing the distance field from a single point source, where average is taken over 100 random samples. Root-mean-square errors are computed with respect to the exact polyhedral solution from [Qin et al. 2016].


# Build
If you have Ninja installed.
```bash
python3 scripts/build.py release
```
Otherwise
```
mkdir build
cd build
cmake ..
make
```