#include "klaplace.h"
#include "vtkio.h"
#include "kgeometry.h"
#include "kvolume.h"

using namespace std;
using namespace pi;


static vtkIO vio;



void processGeometryOptions(Options& opts) {
    opts.addOption("-conv", "File type conversion", SO_NONE);
}


void processGeometryCommands(Options& opts, StringVector& args) {
    if (opts.GetBool("-conv")) {
        vtkDataSet* ds = vio.readDataFile(args[0]);
        vio.writeFile(args[1], ds);
    }
}

int main(int argc, char* argv[]) {
	clock_t t1 = clock();

    Options opts;
    opts.addOption("-h", "print help message", SO_NONE);

	processGeometryOptions(opts);
    processVolumeOptions(opts);

    StringVector args = opts.ParseOptions(argc, argv, NULL);

    if (argc == 1 || opts.GetBool("-h")) {
        cout << "## *kmesh* Usage" << endl;
        opts.PrintUsage();
        return 0;
    }

    processVolumeCommands(opts, args);
	processGeometryCommands(opts, args);
	
	clock_t t2 = clock();
	cout << "Elapsed Time: " << (t2-t1)*(1e-3) << " ms" << endl;
	
    return 0;
}
