# klaplace: Laplacian surface correspondence using an RK45 solver

## Description
The cortical thickness of the mammalian brain is an important morphological characteristic that can be used to investigate and observe the brain’s developmental changes that might be caused by biologically toxic substances such as ethanol or cocaine. Although various cortical thickness analysis methods have been proposed that are applicable for human brain and have developed into well-validated open-source software packages, cortical thickness analysis methods for rodent brains have not yet become as robust and accurate as those designed for human brains. Based on a previously proposed cortical thickness measurement pipeline for rodent brain analysis,1 we present an enhanced cortical thickness pipeline in terms of accuracy and anatomical consistency. First, we propose a Lagrangian-based computational approach in the thickness measurement step in order to minimize local truncation error using the fourth-order Runge-Kutta method. Second, by constructing a line object for each streamline of the thickness measurement, we can visualize the way the thickness is measured and achieve sub-voxel accuracy by performing geometric post-processing. Last, with emphasis on the importance of an anatomically consistent partial differential equation (PDE) boundary map, we propose an automatic PDE boundary map generation algorithm that is specific to rodent brain anatomy, which does not require manual labeling. The results show that the proposed cortical thickness pipeline can produce statistically significant regions that are not observed in the previous cortical thickness analysis pipeline.

## Requirements for build
* VTK5 or above (<a href="https://vtk.org/download/">download</a>)

## Usage
* Input: two surface models (`vtp` or `vtk`)
* Output: deformed source surface to match target surface. The resulting surface will have the same number of vertices as that of the target surface
To find a shape correspondence between `source` and `target`:
```
$ klaplace [-dims grid_spacing] <source> <target> -surfaceCorrespondence <output>
```
The postfix of `_warpedMesh.vtp` will be affixed to `output`. To convert `.vtp` into `.vtk`:
```
$ klaplace -conv output_warpedMesh.vtp output_warpedMesh.vtk
```

## References
* Lee, Joohwi, Sun Hyung Kim, Ipek Oguz, and Martin Styner. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4825173/">"Enhanced cortical thickness measurements for rodent brains via lagrangian-based rk4 streamline computation."</a> In Medical Imaging 2016: Image Processing, vol. 9784, p. 97840B. International Society for Optics and Photonics, 2016.</li>

## Acknowledgment
* Credit to Joohwi Lee for all his implementation of this tool.
* This tool is a modified version of klaplace specifically desinged for <a href="https://github.com/ilwoolyu/LocalGyrificationIndex">Local Gyrification Index</a>. See `master` branch for the original implementation.
