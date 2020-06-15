//
//  kvolume.cpp
//  ktools
//
//  Created by Joowhi Lee on 9/3/15.
//
//

#include "kvolume.h"
#include "kgeometry.h"
#include "vtkio.h"
#include "piOptions.h"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkGradientFilter.h>
#include <vtkMath.h>
#include <vtkThresholdPoints.h>
#include <vtkCleanPolyData.h>
#include <vtkModifiedBSPTree.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkGenericCell.h>
#include <vtkPolyDataNormals.h>
#include <vtkImageData.h>
#include <vtkImageStencil.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageToStructuredGrid.h>
#include <vtkStreamTracer.h>

#include <ctime>
#include <algorithm>
#include <omp.h>

using namespace pi;
using namespace std;
using namespace std::tr1;

static vtkIO vio;



void findNeighborPoints(vtkCell* cell, vtkIdType pid, set<vtkIdType>& nbrs) {
    for (size_t j = 0; j < cell->GetNumberOfPoints(); j++) {
        vtkIdType cellPt = cell->GetPointId(j);
        if (pid != cellPt) {
            nbrs.insert(cellPt);
        }
    }
}

vtkDataSet* createGrid(vtkPolyData* osurf, vtkPolyData* isurf, const int dims, size_t& insideCountOut) {
	
	// compute the common voxel space
	struct GridCreate {
		int dim[3];
		double center[3];
		double spacing[3];
		int extent[6];
		double origin[3];
		
		GridCreate(double* bounds, const int dims) {
			double maxbound = max(bounds[1]-bounds[0], max(bounds[3]-bounds[2], bounds[5]-bounds[4]));
			
			center[0] = (bounds[1]+bounds[0])/2.0;
			center[1] = (bounds[3]+bounds[2])/2.0;
			center[2] = (bounds[5]+bounds[4])/2.0;
			
			double gridSpacing = maxbound / dims;
			spacing[0] = spacing[1] = spacing[2] = gridSpacing;
			
			dim[0] = (bounds[1]-bounds[0])/gridSpacing;
			dim[1] = (bounds[3]-bounds[2])/gridSpacing;
			dim[2] = (bounds[5]-bounds[4])/gridSpacing;
			
			extent[0] = extent[2] = extent[4] = 0;
			extent[1] = dim[0] + 6;
			extent[3] = dim[1] + 6;
			extent[5] = dim[2] + 6;
			
			origin[0] = bounds[0] + gridSpacing*-3;
			origin[1] = bounds[2] + gridSpacing*-3;
			origin[2] = bounds[4] + gridSpacing*-3;
			
			cout << "Grid Dimension: " << dims << "; Grid Spacing: " << gridSpacing << endl;
		}
		
		
		void createImage(vtkImageData* im) {
			im->SetSpacing(spacing);
			im->SetExtent(extent);
			im->SetOrigin(origin);
			im->AllocateScalars(VTK_INT, 1);
			im->GetPointData()->GetScalars()->FillComponent(0, 255);
		}
		
		vtkStructuredGrid* createStencil(vtkImageData* im, vtkPolyData* surf) {
			createImage(im);
			
			vtkNew<vtkPolyDataToImageStencil> psten;
			psten->SetInputData(surf);
			psten->SetOutputOrigin(origin);
			psten->SetOutputSpacing(spacing);
			psten->SetOutputWholeExtent(extent);
			psten->Update();
			
			vtkNew<vtkImageStencil> isten;
			isten->SetInputData(im);
			isten->SetStencilData(psten->GetOutput());
			isten->ReverseStencilOff();
			isten->SetBackgroundValue(0);
			isten->Update();
			
			vtkImageData* imgGrid = isten->GetOutput();
			imgGrid->GetPointData()->GetScalars()->SetName("SampledValue");
			
			vtkNew<vtkImageToStructuredGrid> imgToGrid;
			imgToGrid->SetInputData(imgGrid);
			imgToGrid->Update();
			
			vtkStructuredGrid* output = imgToGrid->GetOutput();
			output->GetPointData()->SetScalars(imgGrid->GetPointData()->GetScalars());
			output->Register(NULL);
			
			return output;
		}
	};
	
	GridCreate gc(osurf->GetBounds(), dims);
	
	vtkNew<vtkImageData> gim;
	vtkNew<vtkImageData> wim;
	
	vtkStructuredGrid* goim = gc.createStencil(gim.GetPointer(), osurf);
	vtkStructuredGrid* woim = gc.createStencil(wim.GetPointer(), isurf);
	
	struct BoundaryCheck {
		size_t subtract(vtkDataSet* aim, vtkDataSet* bim) {
			if (aim->GetNumberOfPoints() != bim->GetNumberOfPoints()) {
				cout << "can't process: the number of points are different!" << endl;
				return 0;
			}
			
			vtkDataArray* aarr = aim->GetPointData()->GetScalars();
			vtkDataArray* barr = bim->GetPointData()->GetScalars();
			
			size_t insideCount = 0;
			for (size_t j = 0; j < aim->GetNumberOfPoints(); j++) {
				int p = aarr->GetTuple1(j);
				int q = barr->GetTuple1(j);
				int o = 700;
				if (p == 255 && q != 255) {
					o = 1;
					insideCount ++;
				} else if (p == 255 && q == 255){
					o = 300;
				}
				aarr->SetTuple1(j, o);
			}
			return insideCount;
		}
		
		void checkSurface(vtkStructuredGrid* grid, vtkPolyData* isurf,  vtkPolyData* osurf) {
			vtkNew<vtkPointLocator> gridLoc;
			gridLoc->SetDataSet(grid);
			gridLoc->BuildLocator();

			vtkIntArray* sampledValue = vtkIntArray::SafeDownCast(grid->GetPointData()->GetScalars());
			
			const size_t nPoints = isurf->GetNumberOfPoints();
			size_t cnt = 0;
			
			for (size_t j = 0; j < nPoints; j++) {
				double p[3];
				
				isurf->GetPoint(j, p);
				vtkIdType pid = gridLoc->FindClosestPoint(p);
				
				int sample = sampledValue->GetValue(pid);
				if (sample == 300) {
					sampledValue->SetValue(pid, 1);
					cnt++;
				}
			}
			cout << "# of inside boundary correction: " << cnt << endl;
			
			cnt = 0;
			const size_t nPoints2 = osurf->GetNumberOfPoints();
			for (size_t j = 0; j < nPoints2; j++) {
				double p[3];
				
				osurf->GetPoint(j, p);
				vtkIdType pid = gridLoc->FindClosestPoint(p);
				
				int sample = sampledValue->GetValue(pid);
				if (sample == 700) {
					sampledValue->SetValue(pid, 1);
					cnt++;
				}
			}
			cout << "# of outside boundary correction: " << cnt << endl;
		}

	};
	
	BoundaryCheck bc;
	insideCountOut = bc.subtract(goim, woim);
	bc.checkSurface(goim, isurf, osurf);
	woim->Delete();
	
	return goim;
}



// create a structured grid with the size of input
// convert the grid to polydata
// create the intersection between the grid and the polydata
vtkDataSet* runFillGrid(Options& opts, StringVector& args) {
	vtkDataSet* grid = NULL;
	if (opts.GetBool("-humanBrain")) {
		string outputFile = args[2];
		
		vtkIO vio;
		vtkPolyData* osurf = vio.readFile(args[0]);
		vtkPolyData* isurf = vio.readFile(args[1]);
		
		size_t insideCountOut = 0;
		grid = createGrid(osurf, isurf, opts.GetStringAsInt("-dims", 100), insideCountOut);
		//vio.writeFile(outputFile, grid);
		
		cout << "Inside Voxels: " << insideCountOut << endl;
	} else {
		cout << "Not supported yet!" << endl;
//		
//		string inputFile = args[0];
//		string outputFile = args[1];
//		
//		vtkIO vio;
//		vtkPolyData* input = vio.readFile(inputFile);
//		int insideCount = 0;
//		vtkDataSet* output = createGridForSphereLikeObject(input, insideCount, 100, false);
//		vio.writeFile(outputFile, output);
//		cout << "Inside Voxels: " << insideCount << endl;
	}
	return grid;
}


// Compute Laplace PDE based on the adjacency list and border
void computeLaplacePDE(vtkDataSet* data, const double low, const double high, const int nIters, const double dt, vtkPolyData* surfaceData = NULL) {
	
	if (data == NULL) {
		cout << "Data input is NULL" << endl;
		return;
	}
	
	class LaplaceGrid {
	public:
		double low;
		double high;
		double dt;
		vtkDataSet* dataSet;
		vtkPolyData* samplePoints;
		
		vector<vtkIdType> solutionDomain;
		vtkIntArray* boundaryCond;
		vtkPolyData* boundarySurface;
		
		vtkDoubleArray* solution;
		vtkDoubleArray* tmpSolution;
		vtkDataArray* laplaceGradient;
		vtkDoubleArray* laplaceGradientNormals;
		
		
		Geometry geom;
		Geometry::NeighborList nbrs;
		
		
		
		LaplaceGrid(double l, double h, double d, vtkDataSet* ds, vtkPolyData* pd= NULL): low(l), high(h), dt(d), dataSet(ds), boundarySurface(pd) {
			cout << "geometry edge extraction ... " << flush;
			geom.extractNeighbors(ds, nbrs);
			cout << " done " << endl;
			
			// check boundary points
			boundaryCond = vtkIntArray::SafeDownCast(ds->GetPointData()->GetArray("SampledValue"));
			if (boundaryCond == NULL) {
				throw runtime_error("No scalar values for BoundaryPoints");
			}
			
			initializeSolution();
		}
		
		void initializeSolution() {
			cout << "initializing solution grid ... " << flush;
			// low-value 2
			// high-value 1
			solution = vtkDoubleArray::New();
			solution->SetName("LaplacianSolution");
			solution->SetNumberOfComponents(1);
			solution->SetNumberOfTuples(boundaryCond->GetNumberOfTuples());
			solution->FillComponent(0, 0);
			
			tmpSolution = vtkDoubleArray::New();
			tmpSolution->SetName("LaplacianSolution");
			tmpSolution->SetNumberOfComponents(1);
			tmpSolution->SetNumberOfTuples(boundaryCond->GetNumberOfTuples());
			tmpSolution->FillComponent(0, 0);
			
			const size_t nPts = boundaryCond->GetNumberOfTuples();
			for (size_t j = 0; j < nPts; j++) {
				int domain = boundaryCond->GetValue(j);
				double uValue = 0;
				if (domain == 700) {
					// high
					uValue = high;
				} else if (domain == 300){
					// low
					uValue = low;
				} else if (domain == 1) {
					uValue = 0;
					solutionDomain.push_back(j);
				}
				solution->SetValue(j, uValue);
				tmpSolution->SetValue(j, uValue);
			}
			cout << "# of points: " << solutionDomain.size() << endl;
		}
		
		void computeStep() {
			const size_t nPts = solutionDomain.size();
			for (size_t j = 0; j < nPts; j++) {
				vtkIdType centerId = solutionDomain[j];
				Geometry::Neighbors& edgeMap = nbrs[centerId];
				Geometry::Neighbors::iterator iter = edgeMap.begin();
				
				double u = 0;
				double nNbrs = 0;
				for (; iter != edgeMap.end(); iter++) {
					const double du = solution->GetValue(*iter);
					u += du;
					nNbrs ++;

					//                    cout << iter->second.axisAligned << endl;
				}
				u = u / nNbrs;
				tmpSolution->SetValue(centerId, u);
			}
			
			vtkDoubleArray* swapTmp = tmpSolution;
			tmpSolution = solution;
			solution = swapTmp;
//			memcpy(solution->WritePointer(0, nPts), tmpSolution->GetVoidPointer(0), sizeof(double) * nTuples);
//			solution->DeepCopy(tmpSolution);
		}
		
		
		void computeNormals(vtkDataSet* data) {
			/*
			 vtkNew<vtkCellDerivatives> deriv;
			 deriv->SetInput(data);
			 deriv->SetVectorModeToComputeGradient();
			 deriv->Update();
			 vtkDataSet* derivOut = deriv->GetOutput();
			 derivOut->GetCellData()->SetActiveVectors("ScalarGradient");
			 vtkDataArray* scalarGradient = deriv->GetOutput()->GetCellData()->GetArray("ScalarGradient");
			 scalarGradient->SetName("LaplacianGradient");
			 */
			
			vtkNew<vtkGradientFilter> gradFilter;
			gradFilter->SetInputData(data);
			gradFilter->SetInputScalars(vtkDataSet::FIELD_ASSOCIATION_POINTS, "LaplacianSolution");
			gradFilter->SetResultArrayName("LaplacianGradient");
			gradFilter->Update();
			laplaceGradient = gradFilter->GetOutput()->GetPointData()->GetArray("LaplacianGradient");
			laplaceGradient->Register(NULL);
			
			laplaceGradientNormals = vtkDoubleArray::New();
			laplaceGradientNormals->SetName("LaplacianGradientNorm");
			laplaceGradientNormals->SetNumberOfComponents(3);
			laplaceGradientNormals->SetNumberOfTuples(laplaceGradient->GetNumberOfTuples());
			
			const size_t nPts = laplaceGradientNormals->GetNumberOfTuples();
			for (size_t j = 0; j < nPts; j++) {
				double* vec = laplaceGradient->GetTuple3(j);
				double norm = vtkMath::Norm(vec);
				if (norm > 1e-10) {
					laplaceGradientNormals->SetTuple3(j, vec[0]/norm, vec[1]/norm, vec[2]/norm);
				} else {
					laplaceGradientNormals->SetTuple3(j, 0, 0, 0);

				}
			}
			
			data->GetPointData()->AddArray(laplaceGradient);
			data->GetPointData()->SetVectors(laplaceGradientNormals);
		}
		
		void computeExteriorNormals(vtkPolyData* boundarySurface, const double radius = .1) {
			if (boundarySurface == NULL) {
				return;
			}
			vtkNew<vtkPolyDataNormals> normalsFilter;
			normalsFilter->SetInputData(boundarySurface);
			normalsFilter->ComputeCellNormalsOn();
			normalsFilter->ComputePointNormalsOn();
			normalsFilter->Update();
			vtkFloatArray* cellNormals = vtkFloatArray::SafeDownCast(normalsFilter->GetOutput()->GetCellData()->GetNormals());
			
			vtkNew<vtkCellLocator> cloc;
			cloc->SetDataSet(boundarySurface);
			cloc->AutomaticOn();
			cloc->BuildLocator();
			
			dataSet->GetPointData()->SetActiveScalars("SampledValue");
			
			vtkNew<vtkThresholdPoints> threshold;
			threshold->SetInputData(dataSet);
			threshold->ThresholdByUpper(250);
			threshold->Update();
			vtkDataSet* inoutBoundary = threshold->GetOutput();
			vtkIntArray* inoutBoundaryCond = vtkIntArray::SafeDownCast(inoutBoundary->GetPointData()->GetArray("SampledValue"));
			
			
			vtkNew<vtkPointLocator> ploc;
			ploc->SetDataSet(inoutBoundary);
			ploc->AutomaticOn();
			ploc->BuildLocator();
			
			const size_t nPts = dataSet->GetNumberOfPoints();
			for (size_t j = 0; j < nPts; j++) {
				int domain = boundaryCond->GetValue(j);
				if (domain == 700 || domain == 300 || domain == 0) {
					double x[3] = { 0, }, closestPoint[3] = { 0, };
					vtkIdType cellId = -1;
					int subId = 0;
					double dist2 = -1;
					dataSet->GetPoint(j, x);
					vtkNew<vtkGenericCell> closestCell;
					cloc->FindClosestPointWithinRadius(x, radius, closestPoint, cellId, subId, dist2);
					
					float cellNormal[3];
#if VTK_MAJOR_VERSION < 7
					cellNormals->GetTupleValue(cellId, cellNormal);
#else
					cellNormals->GetTypedTuple(cellId, cellNormal);
#endif
					cellNormal[0] = 0;
					vtkMath::Normalize(cellNormal);
					
					if (domain == 0) {
						vtkIdType xId = ploc->FindClosestPoint(x);
						domain = inoutBoundaryCond->GetValue(xId);
						assert(domain == 300 || domain == 700);
					}
					
					
					if (domain == 300) {
						laplaceGradientNormals->SetTuple3(j, -cellNormal[0], -cellNormal[1], -cellNormal[2]);
					} else {
						laplaceGradientNormals->SetTuple3(j, cellNormal[0], cellNormal[1], cellNormal[2]);
					}
				}
			}
		}
	};
	
	
	LaplaceGrid grid(low, high, dt, data, surfaceData);
	
	clock_t t1 = clock();
	
	// main iteration loop
	for (size_t i = 1; i <= nIters; i++) {
		/*if (i%500 == 0) {
			cout << "iteration: " << i << "\t";
			clock_t t2 = clock();
			cout << (double) (t2-t1) / CLOCKS_PER_SEC * 1000 << " ms;" << endl;
			t1 = t2;
		}*/
		grid.computeStep();
	}
	clock_t t2 = clock();
	cout << (double) (t2-t1) / CLOCKS_PER_SEC << " sec;" << endl;
	
	
	// return the solution
	data->GetPointData()->AddArray(grid.solution);
	grid.computeNormals(data);
//	grid.computeExteriorNormals(surfaceData);
}



/// @brief perform a line clipping to fit within the object
bool performLineClipping(vtkPolyData* streamLines, vtkModifiedBSPTree* tree, int lineId, vtkCell* lineToClip, vtkPoints* outputPoints, vtkCellArray* outputLines, double &length) {
	
	/// - Iterate over all points in a line
	vtkIdList* ids = lineToClip->GetPointIds();
	
	/// - Identify a line segment included in the line
	bool foundEndpoint = false;
	std::vector<vtkIdType> idList;

	double x[3] = {-1,-1,-1};
	int j = ids->GetNumberOfIds() - 1;
	for (; j >= 2 && !foundEndpoint; j--) {    // IL: inverse search for fast clipping
		double p1[3], p2[3];
		streamLines->GetPoint(ids->GetId(j-1), p1);
		streamLines->GetPoint(ids->GetId(j), p2);
		
		int subId;
		double t = 0;
		
		double pcoords[3] = { -1, };
		foundEndpoint = tree->IntersectWithLine(p1, p2, 0.01, t, x, pcoords, subId);
		//cout << testLine << "; " << x[0] << "," << x[1] << "," << x[2] << endl;
	}
	
	if (ids->GetNumberOfIds() > 2) {
		j++;
		if (!foundEndpoint) {
			streamLines->GetPoint(ids->GetId(j), x);
		}
		int b = max(0,j-2);
		double p1[3], p2[3];
		streamLines->GetPoint(ids->GetId(b), p1);
		idList.push_back(outputPoints->GetNumberOfPoints());
		outputPoints->InsertNextPoint(p1);
		for (int k = b+1; k < j; k++) {
			streamLines->GetPoint(ids->GetId(k), p2);
			idList.push_back(outputPoints->GetNumberOfPoints());
			outputPoints->InsertNextPoint(p2);
			//length += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
			p1[0] = p2[0]; p1[1] = p2[1]; p1[2] = p2[2];
		}
		idList.push_back(outputPoints->GetNumberOfPoints());
		outputPoints->InsertNextPoint(x);
		//length += sqrt(vtkMath::Distance2BetweenPoints(p1, x));
	}
	outputLines->InsertNextCell(idList.size(), &idList[0]);

	return foundEndpoint;
}

bool performLineClippingSimple(vtkPolyData* streamLines, vtkModifiedBSPTree* tree, vtkCell* lineToClip, vtkPoints* outputPoints, int &numPoints) {
	/// - Iterate over all points in a line
	vtkIdList* ids = lineToClip->GetPointIds();

	/// - Identify a line segment included in the line
	bool foundEndpoint = false;

	double x[3] = {-1,-1,-1};
	int j = ids->GetNumberOfIds() - 1;
	for (; j >= 2 && !foundEndpoint; j--) {    // IL: inverse search for fast clipping
		double p1[3], p2[3];
		streamLines->GetPoint(ids->GetId(j-1), p1);
		streamLines->GetPoint(ids->GetId(j), p2);

		int subId;
		double t = 0;

		double pcoords[3] = { -1, };
		foundEndpoint = tree->IntersectWithLine(p1, p2, 0.01, t, x, pcoords, subId);
		//cout << testLine << "; " << x[0] << "," << x[1] << "," << x[2] << endl;
	}
	numPoints = 0;
	if (ids->GetNumberOfIds() > 2) {
		j++;
		if (!foundEndpoint) {
			streamLines->GetPoint(ids->GetId(j), x);
		}
		int b = max(0,j-2);
		double p[3];
		streamLines->GetPoint(ids->GetId(b), p);
		outputPoints->InsertNextPoint(p);
		numPoints++;
		for (int k = b+1; k < j; k++) {
			streamLines->GetPoint(ids->GetId(k), p);
			outputPoints->InsertNextPoint(p);
			numPoints++;
		}
		outputPoints->InsertNextPoint(x);
		numPoints++;
	}

	return foundEndpoint;
}

vtkPolyData* performStreamTracerPostProcessing(vtkPolyData* streamLines, vtkPolyData* seedPoints, vtkPolyData* destinationSurface) {
	
	const size_t nInputPoints = seedPoints->GetNumberOfPoints();
	
	// remove useless pointdata information
	streamLines->GetPointData()->Reset();
	streamLines->BuildCells();
	streamLines->BuildLinks();
	
	
	// loop over the cell and compute the length
	int nCells = streamLines->GetNumberOfCells();
	
	/// - Prepare the output as a scalar array
	//    vtkDataArray* streamLineLength = streamLines->GetCellData()->GetScalars("Length");
	
	/// - Prepare the output for the input points
	vtkDoubleArray* streamLineLengthPerPoint = vtkDoubleArray::New();
	streamLineLengthPerPoint->SetNumberOfTuples(nInputPoints);
	streamLineLengthPerPoint->SetName("Length");
	streamLineLengthPerPoint->SetNumberOfComponents(1);
	streamLineLengthPerPoint->FillComponent(0, 0);
	
	vtkIntArray* lineCorrect = vtkIntArray::New();
	lineCorrect->SetName("LineOK");
	lineCorrect->SetNumberOfValues(nInputPoints);
	lineCorrect->FillComponent(0, 0);
	
	seedPoints->GetPointData()->SetScalars(streamLineLengthPerPoint);
	seedPoints->GetPointData()->AddArray(lineCorrect);
	
	cout << "Assigning length to each source vertex ..." << endl;
	vtkDataArray* seedIds = streamLines->GetCellData()->GetScalars("SeedIds");
	if (seedIds) {
		// line clipping
		vtkPoints* outputPoints = vtkPoints::New();
		vtkCellArray* outputCells = vtkCellArray::New();
		
		/// construct a tree locator
		vtkModifiedBSPTree* tree = vtkModifiedBSPTree::New();
		tree->SetDataSet(destinationSurface);
		tree->BuildLocator();
		
		vtkDoubleArray* lengthArray = vtkDoubleArray::New();
		lengthArray->SetName("Length");
		
		vtkIntArray* pointIds = vtkIntArray::New();
		pointIds->SetName("PointIds");
		
		cout << "# of cells: " << nCells << endl;

		int noLines = 0;
		for (int i = 0; i < nCells; i++) {
			int pid = seedIds->GetTuple1(i);
			double length = 0;
			if (pid > -1) {
				vtkCell* line = streamLines->GetCell(i);
				/// - Assume that a line starts from a point on the input mesh and must meet at the opposite surface of the starting point.
				bool lineAdded = performLineClipping(streamLines, tree, i, line, outputPoints, outputCells, length);
				
				if (lineAdded) {
					pointIds->InsertNextValue(pid);
					lengthArray->InsertNextValue(length);
					streamLineLengthPerPoint->SetValue(pid, length);
					lineCorrect->SetValue(pid, 1);
				} else {
					pointIds->InsertNextValue(pid);
					lengthArray->InsertNextValue(length);
					streamLineLengthPerPoint->SetValue(pid, length);
					lineCorrect->SetValue(pid, 2);
					noLines++;
				}
			}
		}
		
		cout << "# of clipping failure: " << noLines << endl;

		streamLines->Delete();

		vtkPolyData* outputStreamLines = vtkPolyData::New();
		outputStreamLines->SetPoints(outputPoints);
		outputStreamLines->SetLines(outputCells);
		outputStreamLines->GetCellData()->AddArray(pointIds);
		outputStreamLines->GetCellData()->AddArray(lengthArray);
		
		vtkCleanPolyData* cleaner = vtkCleanPolyData::New();
		cleaner->SetInputData(outputStreamLines);
		cleaner->Update();
		
		return cleaner->GetOutput();
	} else {
		cout << "Can't find SeedIds" << endl;
		return NULL;
	}
}


vtkPolyData* performStreamTracer(Options& opts, vtkDataSet* inputData, vtkPolyData* inputSeedPoints, vtkPolyData* destSurf, bool zRotate = false) {
    if (inputData == NULL || inputSeedPoints == NULL) {
        cout << "input vector field or seed points is null!" << endl;
        return NULL;
    }
    
    if (destSurf == NULL) {
        cout << "trace destination surface is null" << endl;
        return NULL;
    }
    
	// set active velocity field
	inputData->GetPointData()->SetActiveVectors("LaplacianGradientNorm");
	
	/// - Converting the input points to the image coordinate
	vtkPoints* points = inputSeedPoints->GetPoints();
	cout << "# of seed points: " << points->GetNumberOfPoints() << endl;
	const int nInputPoints = inputSeedPoints->GetNumberOfPoints();
	if (zRotate) {
		for (int i = 0; i < nInputPoints; i++) {
			double p[3];
			points->GetPoint(i, p);
			// FixMe: Do not use a specific scaling factor
			if (zRotate) {
				p[0] = -p[0];
				p[1] = -p[1];
				p[2] = p[2];
			}
			points->SetPoint(i, p);
		}
		inputSeedPoints->SetPoints(points);
	}
	
	
	/// StreamTracer should have a point-wise gradient field
	/// - Set up tracer (Use RK45, both direction, initial step 0.05, maximum propagation 500
	vtkStreamTracer* tracer = vtkStreamTracer::New();
	tracer->SetInputData(inputData);
	tracer->SetSourceData(inputSeedPoints);
	tracer->SetComputeVorticity(false);

	string traceDirection = opts.GetString("-traceDirection", "forward");
	if (traceDirection == "both") {
		tracer->SetIntegrationDirectionToBoth();
	} else if (traceDirection == "backward") {
		tracer->SetIntegrationDirectionToBackward();
		cout << "Backward Integration" << endl;
	} else if (traceDirection == "forward") {
		tracer->SetIntegrationDirectionToForward();
		cout << "Forward Integration" << endl;
	}

//    tracer->SetInterpolatorTypeToDataSetPointLocator();
	tracer->SetInterpolatorTypeToCellLocator();
	tracer->SetMaximumPropagation(5000);
	tracer->SetInitialIntegrationStep(0.01);
//    tracer->SetMaximumIntegrationStep(0.1);
    tracer->SetIntegratorTypeToRungeKutta45();
//    tracer->SetIntegratorTypeToRungeKutta2();

    cout << "Integration Direction: " << tracer->GetIntegrationDirection() << endl;
    cout << "Initial Integration Step: " << tracer->GetInitialIntegrationStep() << endl;
    cout << "Maximum Integration Step: " << tracer->GetMaximumIntegrationStep() << endl;
    cout << "Minimum Integration Step: " << tracer->GetMinimumIntegrationStep() << endl;
    cout << "Maximum Error: " << tracer->GetMaximumError() << endl;
    cout << "IntegratorType: " << tracer->GetIntegratorType() << endl;

    
    tracer->Update();

	vtkPolyData* streamLines = tracer->GetOutput();
//	streamLines->Print(cout);
	
//	vio.writeFile("streamlines.vtp", streamLines);
	
	return streamLines;
}

vtkPolyData* performStreamTracerBatch(Options& opts, vtkDataSet* inputData, vtkPolyData* inputSeedPoints, vtkPolyData* destSurf) {
    if (inputData == NULL || inputSeedPoints == NULL) {
        cout << "input vector field or seed points is null!" << endl;
        return NULL;
    }

    if (destSurf == NULL) {
        cout << "trace destination surface is null" << endl;
        return NULL;
    }

	// set active velocity field
	inputData->GetPointData()->SetActiveVectors("LaplacianGradientNorm");

	/// - Converting the input points to the image coordinate
	vtkPoints* points = inputSeedPoints->GetPoints();
	const int nInputPoints = inputSeedPoints->GetNumberOfPoints();
	cout << "# of seed points: " << nInputPoints << endl;

	// line clipping
	vtkPoints* outputPoints = vtkPoints::New();
	vtkCellArray* outputCells = vtkCellArray::New();

	vtkIntArray* pointIds = vtkIntArray::New();
	pointIds->SetName("PointIds");

	int thread = opts.GetStringAsInt("-thread", 0);
	if (thread <= 0) {
		const char *env = getenv("OMP_NUM_THREADS");
		thread = (env != NULL) ? max(atoi(env), 1) : 1;
	}
	omp_set_num_threads(thread);
	cout << "OpenMP threads: " << thread << endl;

	vtkStreamTracer** tracer_array = new vtkStreamTracer*[thread];
	vtkPolyData** poly_array = new vtkPolyData*[thread];
	vtkPoints** point_array = new vtkPoints*[thread];
	vtkPolyData** streamLines_array = new vtkPolyData*[thread];
	vtkPoints** outputPoints_array = new vtkPoints*[thread];
	vtkIntArray** pointIds_array = new vtkIntArray*[thread];
	vtkIntArray** pointNum_array = new vtkIntArray*[thread];
	vtkModifiedBSPTree** tree_array = new vtkModifiedBSPTree*[thread];

	string traceDirection = opts.GetString("-traceDirection", "forward");
	for (int i = 0; i < thread; i++)
	{
		tracer_array[i] = vtkStreamTracer::New();
		poly_array[i] = vtkPolyData::New();
		point_array[i] = vtkPoints::New();
		outputPoints_array[i] = vtkPoints::New();
		pointIds_array[i] = vtkIntArray::New();
		pointNum_array[i] = vtkIntArray::New();
		// construct a tree locator
		tree_array[i] = vtkModifiedBSPTree::New();
		tree_array[i]->SetDataSet(destSurf);
		tree_array[i]->BuildLocator();

		vtkStreamTracer* tracer = tracer_array[i];
		// tracer->SetInterpolatorTypeToDataSetPointLocator();
		tracer->SetInterpolatorTypeToCellLocator();
		tracer->SetMaximumPropagation(5000);
		tracer->SetInitialIntegrationStep(0.01);
		// tracer->SetMaximumIntegrationStep(0.1);
		tracer->SetIntegratorTypeToRungeKutta45();
		// tracer->SetIntegratorTypeToRungeKutta2();
		tracer->SetInputData(inputData);
		tracer->SetComputeVorticity(false);

		if (traceDirection == "both") {
			tracer->SetIntegrationDirectionToBoth();
		} else if (traceDirection == "backward") {
			tracer->SetIntegrationDirectionToBackward();
			// cout << "Backward Integration" << endl;
		} else if (traceDirection == "forward") {
			tracer->SetIntegrationDirectionToForward();
			// cout << "Forward Integration" << endl;
		}

		vtkPolyData* poly = poly_array[i];
		tracer->SetSourceData(poly);
		vtkPoints* point = point_array[i];
		poly->SetPoints(point);
		streamLines_array[i] = tracer->GetOutput();
	}

	int noLines = 0;
	// int partition = 8;
	// int batch = (nInputPoints / partition + (nInputPoints % partition > 0));
	int batch = opts.GetStringAsInt("-batchSeed", 20000);
	if (batch < 1) batch = 20000;

	for (int i = 0; i < nInputPoints; i += batch) {
		cout << "Batch " << (i / batch + 1) << " ... ";
		fflush(stdout);
		/// StreamTracer should have a point-wise gradient field
		/// - Set up tracer (Use RK45, both direction, initial step 0.05, maximum propagation 500
		int size = min(batch, nInputPoints - i);
		int min_batch = (size / thread + (size % thread > 0));
		#pragma omp parallel for
		for (int n = 0; n < thread; n++)
		{
			vtkStreamTracer* tracer = tracer_array[n];
			vtkPoints* point = point_array[n];
			vtkPolyData* streamLines = streamLines_array[n];
			point->Reset();
			pointIds_array[n]->Reset();
			outputPoints_array[n]->Reset();
			pointNum_array[n]->Reset();
			for (int j = min_batch * n; j < min(min_batch * (n + 1), size); j++) {
				double p[3];
				points->GetPoint(i + j, p);
				point->InsertNextPoint(p);
			}
			tracer->Update();

			streamLines->GetPointData()->Reset();
			streamLines->BuildCells();
			streamLines->BuildLinks();

			vtkDataArray* seedIds = streamLines->GetCellData()->GetScalars("SeedIds");
			if (seedIds) {
				int nCells = streamLines->GetNumberOfCells();
				for (int j = 0; j < nCells; j++) {
					int pid = seedIds->GetTuple1(j);
					double length = 0;
					if (pid > -1) {
						pid += (i + n * min_batch);
						vtkCell* line = streamLines->GetCell(j);
						int numPoints;
						/// - Assume that a line starts from a point on the input mesh and must meet at the opposite surface of the starting point.
						bool lineAdded = performLineClippingSimple(streamLines, tree_array[n], line, outputPoints_array[n], numPoints);
						pointIds_array[n]->InsertNextValue(pid);
						pointNum_array[n]->InsertNextValue(numPoints);
						if (!lineAdded) {
							#pragma omp atomic
							noLines++;
						}
					}
				}
			} else {
				cout << "Can't find SeedIds" << endl;
				exit(1);
			}
		}
		for (int n = 0; n < thread; n++) {
			int id = 0;
			std::vector<vtkIdType> idList;
			for (int j = 0; j < outputPoints_array[n]->GetNumberOfPoints(); j++) {
				idList.push_back(outputPoints->GetNumberOfPoints());
				if (idList.size() == pointNum_array[n]->GetValue(id)) {
					outputCells->InsertNextCell(idList.size(), &idList[0]);
					idList.clear();
					id++;
				}
				outputPoints->InsertNextPoint(outputPoints_array[n]->GetPoint(j));
			}
			for (int j = 0; j < pointIds_array[n]->GetNumberOfValues(); j++) {
				pointIds->InsertNextValue(pointIds_array[n]->GetValue(j));
			}
		}
		cout << "done" << endl;
	}
	cout << "# of clipping failure: " << noLines << endl;
	vtkPolyData* outputStreamLines = vtkPolyData::New();
	outputStreamLines->SetPoints(outputPoints);
	outputStreamLines->SetLines(outputCells);
	outputStreamLines->GetCellData()->AddArray(pointIds);

	vtkCleanPolyData* cleaner = vtkCleanPolyData::New();
	cleaner->SetInputData(outputStreamLines);
	cleaner->Update();

	return cleaner->GetOutput();
}

void runPrintTraceCorrespondence_(Options& opts, string inputMeshName, vtkDataSet* strmesh, string outputWarpedMeshName, vtkPolyData* srcmesh) {
	vtkIO vio;
	
	if (srcmesh == NULL) {
		srcmesh = vio.readFile(inputMeshName);
	}
	
	vtkNew<vtkPolyData> warpedMesh;
	warpedMesh->DeepCopy(srcmesh);
	
	srcmesh->ComputeBounds();
	
	double center[3];
	srcmesh->GetCenter(center);
	
	int traceDirection = vtkStreamTracer::FORWARD;
	string dir = opts.GetString("-traceDirection", "forward");
	if (dir == "backward") {
		traceDirection = vtkStreamTracer::BACKWARD;
	} else if (dir == "both") {
		traceDirection = vtkStreamTracer::BOTH;
	}
	
	vtkNew<vtkDoubleArray> pointArr;
	pointArr->SetName("SourcePoints");
	pointArr->SetNumberOfComponents(3);
	pointArr->SetNumberOfTuples(srcmesh->GetNumberOfPoints());
	
	vtkNew<vtkDoubleArray> destPointArr;
	destPointArr->SetName("DestinationPoints");
	destPointArr->SetNumberOfComponents(3);
	destPointArr->SetNumberOfTuples(srcmesh->GetNumberOfPoints());
	
	vtkNew<vtkPoints> warpedPoints;
	warpedPoints->DeepCopy(srcmesh->GetPoints());
	
	vtkNew<vtkPointLocator> ploc;
	ploc->SetDataSet(srcmesh);
	ploc->SetTolerance(0);
	ploc->BuildLocator();

	vtkDataArray* seedIds = strmesh->GetCellData()->GetArray("PointIds");
	
	const size_t nCells = strmesh->GetNumberOfCells();
	for (size_t j = 0; j < nCells; j++) {
		vtkCell* cell = strmesh->GetCell(j);
		const size_t nPts = cell->GetNumberOfPoints();
		if (nPts < 2) {
			continue;
		}
		vtkIdType e = cell->GetPointId(nPts-1);
		
		double qe[3];
		strmesh->GetPoint(e, qe);
		
		vtkIdType seedId = (vtkIdType) seedIds->GetTuple1(j);
		warpedPoints->SetPoint(seedId, qe);
	}

	struct InterpolateBrokenPoints {
		InterpolateBrokenPoints(vtkPolyData* surf, vtkPoints* warpedPoints, vtkDataArray* seedIds) {
			// identify broken points
			vector<vtkIdType> brokenPoints;
			vtkIdType z = 0;
			for (size_t j = 0; j < seedIds->GetNumberOfTuples(); j++,z++) {
				vtkIdType y = seedIds->GetTuple1(j);
				while (z < y) {
					brokenPoints.push_back(z++);
				}
			}
			
			// find neighbors and compute interpolatead points
			vtkNew<vtkIdList> cellIds;
			set<vtkIdType> nbrs;
			for (size_t j = 0; j < brokenPoints.size(); j++) {
				vtkIdType pid = brokenPoints[j];
				cellIds->Reset();
				surf->GetPointCells(pid, cellIds.GetPointer());
				nbrs.clear();
				// find neighbor points
				for (size_t k = 0; k < cellIds->GetNumberOfIds(); k++) {
					vtkCell* cell = surf->GetCell(k);
					findNeighborPoints(cell, pid, nbrs);
				}
				// average neighbor points
				double p[3] = {0,}, q[3] = {0,};
				set<vtkIdType>::iterator it = nbrs.begin();
				for (; it != nbrs.end(); it++) {
					//if (find(brokenPoints.begin(), brokenPoints.end(), *it) == brokenPoints.end()) {
					if (!binary_search(brokenPoints.begin(), brokenPoints.end(), *it)) {	// IL: binary search
						warpedPoints->GetPoint(*it, q);
						vtkMath::Add(p, q, p);
					} else {
						cout << "broken neighbor!! " << pid << "," << *it << endl;
					}
				}
				p[0]/=nbrs.size();
				p[1]/=nbrs.size();
				p[2]/=nbrs.size();
				warpedPoints->SetPoint(pid, p);
			}
		}

	};
	
	warpedMesh->SetPoints(warpedPoints.GetPointer());
	InterpolateBrokenPoints(warpedMesh.GetPointer(), warpedPoints.GetPointer(), seedIds);
	
	//	warpedMesh->Print(cout);
//	warpedMesh->GetPointData()->SetVectors(pointArr.GetPointer());
//	warpedMesh->GetPointData()->AddArray(sphrCoord.GetPointer());
	vio.writeFile(outputWarpedMeshName, warpedMesh.GetPointer());
	
//	if (args.size() > 3) {
//		srcmesh->GetPointData()->SetVectors(destPointArr.GetPointer());
//		srcmesh->GetPointData()->AddArray(sphrCoord.GetPointer());
//		vio.writeFile(args[3], srcmesh);
//	}
}

void runPrintTraceCorrespondence(Options& opts, string inputMeshName, string inputStreamName, string outputWarpedMeshName, vtkPolyData* srcmesh = NULL) {
	vtkDataSet* strmesh = vio.readDataFile(inputStreamName);
	runPrintTraceCorrespondence_(opts, inputMeshName, strmesh, outputWarpedMeshName, srcmesh);
}

void runSurfaceCorrespondence(Options& opts, StringVector& args) {
	// runEnclosingSphere
	string inputObj1 = args[0];
    string inputObj2 = args[1];
    string prefix = args[2];

    if (inputObj1 == "" || inputObj2 == "") {
        cout << "-surfaceCorrespondence option needs two inputs" << endl;
    }
    if (prefix == "") {
        prefix = "surface_correspondence";
    }

	string outputGrid = prefix + "_grid.vts";
	string outputField = prefix + "_field.vts";
	string outputStream = prefix + "_stream.vtp";
	string outputMesh = prefix + "_warpedMesh.vtp";
	string outputObj = prefix + "_object.vtp";

    //cout << "Output grid: " << outputGrid << endl;
    //cout << "Output laplacian field: " << outputField << endl;
    //cout << "Output streamlines: " << outputStream << endl;
    cout << "Output warped mesh: " << outputMesh << endl;
	
    string inputFieldFile = opts.GetString("-inputField");
    
    
    vtkDataSet* laplaceField = NULL;
    
    if (opts.HasString("-inputField")) {
        cout << "Reading " << inputFieldFile << flush;
        laplaceField = vio.readDataFile(inputFieldFile);
        cout << " done." << endl;
    } else {
        // create uniform grid for a FDM model
        StringVector fillGridArgs;
        fillGridArgs.push_back(inputObj2);
        fillGridArgs.push_back(inputObj1);
        fillGridArgs.push_back(outputGrid);
        opts.SetBool("-humanBrain", true);
        laplaceField = runFillGrid(opts, fillGridArgs);
        
        // compute laplace map
        //laplaceField = vio.readDataFile(outputGrid);
        
        const double dt = opts.GetStringAsReal("-dt", 0.125);
        const int numIter = opts.GetStringAsInt("-iter", 10000);
        
        computeLaplacePDE(laplaceField, 0, 10000, numIter, dt);
        //vio.writeFile(outputField, laplaceField);
    }
	
	
	
	vtkIO vio;
	vtkPolyData* inputData = vio.readFile(inputObj1);
    vtkPolyData* inputData2 = vio.readFile(inputObj2);

    if (inputData == NULL) {
        cout << inputObj1 << " is null" << endl;
        return;
    }
    if (inputData2 == NULL) {
        cout << inputObj2 << " is null" << endl;
        return;
    }

	if (opts.GetString("-traceDirection") == "backward") {
		/*vtkPolyData* streams = performStreamTracer(opts, laplaceField, inputData2, inputData);
		laplaceField->Delete();
		streams = performStreamTracerPostProcessing(streams, inputData2, inputData);*/
		vtkPolyData* streams = performStreamTracerBatch(opts, laplaceField, inputData2, inputData);
		laplaceField->Delete();
		/*vio.writeFile(outputStream, streams);
		runPrintTraceCorrespondence(opts, inputObj2, outputStream, outputMesh, inputData2);*/
		runPrintTraceCorrespondence_(opts, inputObj2, streams, outputMesh, inputData2);
	} else {
		/*vtkPolyData* streams = performStreamTracer(opts, laplaceField, inputData, inputData2);
		laplaceField->Delete();
		streams = performStreamTracerPostProcessing(streams, inputData, inputData2);*/
		vtkPolyData* streams = performStreamTracerBatch(opts, laplaceField, inputData, inputData2);
		laplaceField->Delete();
		/*vio.writeFile(outputStream, streams);
		runPrintTraceCorrespondence(opts, inputObj1, outputStream, outputMesh, inputData);*/
		runPrintTraceCorrespondence_(opts, inputObj1, streams, outputMesh, inputData);
	}
	

	//vio.writeFile(outputObj, inputData);
}



void processVolumeOptions(Options& opts) {
	opts.addOption("-dims", "x-y-z dimensions", "-dims 100", SO_REQ_SEP);
	
	opts.addOption("-surfaceCorrespondence", "Construct a surface correspondence between two objects; prefix is used for temporary files", "-surfaceCorrespondence source.vtp destination.vtp prefix", SO_NONE);

	opts.addOption("-batchSeed", "# of seeds to be traced once (default: 20000)", "", SO_REQ_SEP);

	opts.addOption("-thread", "# of OpenMP threads", "", SO_REQ_SEP);
}

void processVolumeCommands(Options& opts, StringVector& args) {
	
	string input1File, outputFile;
	
    if (opts.GetBool("-surfaceCorrespondence")) {
        runSurfaceCorrespondence(opts, args);
    }
}
