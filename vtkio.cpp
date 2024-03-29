//
//  vtkio.cpp
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#include "vtkio.h"
#include <exception>
#include <stdexcept>

#include <vtkNew.h>
#include <vtkDataSet.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkFieldData.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkMNIObjectReader.h>
#include <vtkMNIObjectWriter.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkMath.h>
#include <vtkPointData.h>

static bool endswith(std::string file, std::string ext) {
    int epos = file.length() - ext.length();
    if (epos < 0) {
        return false;
    }
    return file.rfind(ext) == epos;
}

void vtkIO::zrotate(vtkPolyData* p) {
    int np = p->GetNumberOfPoints();
    vtkPoints* points = p->GetPoints();
    for (int i = 0; i < np; i++) {
        double x[3];
        points->GetPoint(i, x);
//        cout << x[0] << "," << x[1] << "," << x[2] << endl;
        x[0] = -x[0];
        x[1] = -x[1];
        points->SetPoint(i, x);
//        cout << x[0] << "," << x[1] << "," << x[2] << endl;
    }
    p->SetPoints(points);
}


/// @brief Read a vtk/vtp file. The file's type is automatically determined by its extension.
vtkPolyData* vtkIO::readFile(std::string file) {
    if (endswith(file, ".vtp")) {
        vtkXMLPolyDataReader* r = vtkXMLPolyDataReader::New();
        r->SetFileName(file.c_str());
        r->Update();
        return r->GetOutput();
    } else if (endswith(file, ".vtk")) {
        vtkPolyDataReader* r = vtkPolyDataReader::New();
        r->SetFileName(file.c_str());
        r->Update();
        return r->GetOutput();
    } else if (endswith(file, ".obj")) {
        vtkMNIObjectReader* r = vtkMNIObjectReader::New();
        r->SetFileName(file.c_str());
        r->Update();
        return r->GetOutput();
	} else {
		throw std::runtime_error("file format not recognized by extesion");
	}
    return NULL;
}


vtkDataArray* vtkIO::findFieldData(vtkPolyData* dataSet, std::string propertyName) {
    vtkDataArray* dataArray = NULL;
    int nArr = dataSet->GetFieldData()->GetNumberOfArrays();
    
    dataSet->Print(cout);

    cout << "number of arrays found: " << nArr << endl;
    for (int j = 0; j < nArr; j++) {
        std::string arrayName(dataSet->GetFieldData()->GetArrayName(j));
        std::cout << "found " << arrayName << std::endl;
        if (propertyName == arrayName) {
            dataArray = dataSet->GetPointData()->GetArray(j);
            break;
        }
    }
    return dataArray;
}


/// @brief Write a vtk/vtp file. The file's type is automatically determined by its extension.
void vtkIO::writeFile(std::string file, vtkDataSet *mesh) {
    if (endswith(file, ".vtp")) {
        vtkXMLPolyDataWriter* w = vtkXMLPolyDataWriter::New();
        w->SetInputData(mesh);
        w->SetFileName(file.c_str());
        w->Write();
        w->Delete();
    } else if (endswith(file, ".vtk")) {
        vtkPolyDataWriter* w = vtkPolyDataWriter::New();
#if VTK_MAJOR_VERSION >= 9
        w->SetFileVersion(42);
#endif
        w->SetInputData(mesh);
        w->SetFileName(file.c_str());
        w->Write();
        w->Delete();
    } else if (endswith(file, ".vtu")) {
        vtkXMLUnstructuredGridWriter* w = vtkXMLUnstructuredGridWriter::New();
        w->SetInputData(mesh);
        w->SetFileName(file.c_str());
        w->SetCompressorTypeToNone();
        w->SetDataModeToAscii();
        w->Write();
        w->Delete();
    } else if (endswith(file, ".vts")) {
        vtkXMLStructuredGridWriter* w = vtkXMLStructuredGridWriter::New();
        w->SetInputData(mesh);
        w->SetFileName(file.c_str());
        w->SetCompressorTypeToNone();
        w->Write();
        w->Delete();
        } else if (endswith(file, ".mhd")) {
            vtkMetaImageWriter* w = vtkMetaImageWriter::New();

            w->SetInputData(vtkImageData::SafeDownCast(mesh));
            w->SetFileName(file.c_str());
            w->Write();
            w->Delete();
        }
        cout << "Write " << file << " done ..." << endl;
}


void vtkIO::writeXMLFile(std::string file, vtkPolyData *mesh) {
    vtkXMLPolyDataWriter* w = vtkXMLPolyDataWriter::New();
    w->SetInputData(mesh);
    w->SetFileName(file.c_str());
    w->SetDataModeToAppended();
    w->EncodeAppendedDataOff();
    w->SetCompressorTypeToZLib();
    w->SetDataModeToBinary();
    w->Write();
    w->Delete();
}
//
//vtkDataArray* vtkIO::getPointArray(vtkDataSet* ds) {
//    if (ds == NULL) {
//        return NULL;
//    }
//
//    vtkDataArray* da = NULL;
//    for (int i = 0; i < ds->GetPointData()->GetNumberOfArrays(); i++) {
//        ds->GetPointData()->GetArray(const char *arrayName)
//    }
//
//    return da;
//}





vtkDataSet* vtkIO::readDataFile(std::string file) {
	if (endswith(file, ".vtp")) {
		vtkXMLPolyDataReader* r = vtkXMLPolyDataReader::New();
		r->SetFileName(file.c_str());
		r->Update();
		return r->GetOutput();
	} else if (endswith(file, ".vts")) {
		vtkXMLStructuredGridReader* r = vtkXMLStructuredGridReader::New();
		r->SetFileName(file.c_str());
		r->Update();
		return r->GetOutput();
	} else if (endswith(file, ".vtu")) {
		vtkXMLUnstructuredGridReader* r = vtkXMLUnstructuredGridReader::New();
		r->SetFileName(file.c_str());
		r->Update();
		return r->GetOutput();
	} else if (endswith(file, ".vtk")) {
		vtkPolyDataReader* r = vtkPolyDataReader::New();
		r->SetFileName(file.c_str());
		r->Update();
		return r->GetOutput();
	} else if (endswith(file, ".obj")) {
		vtkMNIObjectReader* r = vtkMNIObjectReader::New();
		r->SetFileName(file.c_str());
		r->Update();
		return r->GetOutput();
	} else if (endswith(file, ".mhd")) {
		vtkNew<vtkMetaImageReader> r;
		r->SetFileName(file.c_str());
		r->Update();
		vtkDataSet* ds = r->GetOutput();
		ds->Register(NULL);
		return ds;
	}
	cout << "Unknown file format: " << file << endl;
	return NULL;
}

