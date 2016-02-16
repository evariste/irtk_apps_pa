/*
 * polydatafiles2csv.cc
 *
 *  Created on: Feb 13, 2015
 *      Author: paulaljabar
 */

#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
#include <string.h>
//#include <sys/stat.h>
//
#include <irtkImage.h>

#include <irtkMatrix.h>
#include <abcdUtils.h>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include <vtkSphereSource.h>
#include <vtkTensorGlyph.h>

void usage(std::string name) {

  cout << "" << endl;
  cout << "" << endl;
  cout << name
      << " [output.vtk] [points_1.vtk] [points_2.vtk] .... [points_N.vtk] <options> "
      << endl;
  cout << "" << endl;
  cout << "Generate Gaussian ellipsoids for the landmark data given in the files " << endl;
  cout << "[points_1.vtk] .... [points_N.vtk]. I.e. each file is expected to have" << endl;
  cout << "the same number of corresponding landmarks. The ellipsoids indicate the" << endl;
  cout << "spatial spread of each landmark." << endl;
  cout << "" << endl;
  cout << "This is the 3-D analogue of using an interval in 1-D of, say,  " << endl;
  cout << "  [mean - 3*SD, mean + 3*SD]" << endl;
  cout << "to define a range where we expect to find nearly all the data. " << endl;
  cout << "" << endl;
  cout << "The axes of each ellipsoid are aligned with the principle modes of" << endl;
  cout << "variation for the corresponding set of landmarks. The radius of each" << endl;
  cout << "ellipsoid along each axis is proportional to the standard deviation " << endl;
  cout << "of the data along that axis (see -radius option)." << endl;
  cout << "" << endl;

  cout << "Options:" << endl;
  cout << "" << endl;
  cout << "-radius [value]  : Radius along each axis of the ellipsoid. The number" << endl;
  cout << "                   given is multiplied by the standard deviation of the " << endl;
  cout << "                   the data along the axis. Setting [value] to the following"  << endl;
  cout << "                   gives ellipsoids that are expected to capture the given"  << endl;
  cout << "                   fractions of the landmark cloud."  << endl;
  cout << "" << endl;
  cout << "                    value=1              19.9%               " << endl;
  cout << "                    value=2              73.9%               " << endl;
  cout << "                    value=3 (default)    97.1%               " << endl;
  cout << "" << endl;
  cout << "See https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf" << endl;
  cout << "" << endl;
  cout << "" << endl;
  /*
   *
   * The normal distribution in dimension 1 has an expected proprortion of the
   * population given by a number of standard deviations. These can be obtained
   * from the chi squared cumulative distribution function. In MatLab this is
   * chi2cdf:
   *
   * >> [chi2cdf(1,1) chi2cdf(4,1) chi2cdf(9,1)]
   0.6827    0.9545    0.9973

   * In a two dimensional Gaussian, the region can be defined by the ellipsoid
   * which is at a fixed Mahalonibis distance from the mean. The M. distance is
   * (x - mu)^T \sigma^{-1} (x - mu) where \sigma is the covariance matrix.
   * The ellipsoid has radii along its principle axes that can be expressed as multiples
   * of the standard deviation along that axis. Using 1, 2 or 3 multiples
   * of the standard deviation, we expect to capture the following proportions
   * of the population within each of the successively larger ellipses.
   *
   >> [chi2cdf(1^2, 2) chi2cdf(2^2, 2) chi2cdf(3^2, 2)]
   0.3935    0.8647    0.9889

   * In 3-D the numbers are as follows:
   *
   >> [chi2cdf(1^2, 3) chi2cdf(2^2, 3) chi2cdf(3^2, 3)]
   0.1987    0.7385    0.9707
   *
   * https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf
   *
   */
  exit(1);
}

char **inputFilenames = NULL;
char *outputFileName = NULL;

int main(int argc, char **argv) {

  bool ok, verbose;
  double p[3];
  double tensorVec[9];
  double radius = 3.0;
  int i, noOfPointSets, j;
  int a, b;

  std::string exeName = argv[0];

  if (argc < 2)
    usage(exeName);

  int MAX_FILES = 1000;

  // Parse file names and values

  outputFileName = argv[1];
  argc--;
  argv++;

  char **inputFilenames = new char *[MAX_FILES];

  noOfPointSets = 0;
  while ((argc > 1) && (argv[1][0] != '-')) {
    inputFilenames[noOfPointSets] = argv[1];

    if (!is_vtkPolyDataFile(inputFilenames[noOfPointSets])) {
      cerr << "Not a vtk polydata file: " << inputFilenames[noOfPointSets]
          << endl;
      usage(exeName);
    }

    noOfPointSets++;
    argc--;
    argv++;
  }

  if (noOfPointSets < 1) {
    cerr << "No files to read." << endl;
    usage(exeName);
  }

  if (noOfPointSets < 2) {
    cerr << "Only one file. Covariance across landmarks meaningless." << endl;
    exit(1);
  }

  verbose = false;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-verbose") == 0)) {
      argc--;
      argv++;
      verbose = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-radius") == 0)) {
      argc--;
      argv++;
      radius = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  if (verbose) {
    cerr << "Processing " << noOfPointSets << " files " << endl;
    for (i = 0; i < noOfPointSets; i++)
      cerr << "   " << inputFilenames[i] << endl;
    if (outputFileName != NULL)
      cerr << "Writing output to " << outputFileName << endl;
  }

  vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();

  reader->SetFileName(inputFilenames[0]);
  reader->Update();

  vtkSmartPointer<vtkPolyData> inputPoints =
      vtkSmartPointer<vtkPolyData>::New();

  inputPoints = reader->GetOutput();

  int nPts = inputPoints->GetNumberOfPoints();

  FILE *output;

  output = fopen(outputFileName, "w");

  double ***data = NULL;
  double **meanXYZ = NULL;

  data = Allocate(data, 3, nPts, noOfPointSets);

  for (i = 0; i < noOfPointSets; i++) {

    reader->SetFileName(inputFilenames[i]);
    reader->Update();

    inputPoints = reader->GetOutput();
    if (inputPoints->GetNumberOfPoints() != nPts) {
      cerr << "Input files have differing numbers of points." << endl;
      exit(1);
    }

    for (j = 0; j < nPts; j++) {
      // Read all points for current set.
      inputPoints->GetPoint(j, p);
      data[i][j][0] = p[0];
      data[i][j][1] = p[1];
      data[i][j][2] = p[2];
    }

  }

  // mean location for each point
  meanXYZ = Allocate(meanXYZ, 3, nPts);

  for (j = 0; j < nPts; j++) {
    meanXYZ[j][0] = 0;
    meanXYZ[j][1] = 0;
    meanXYZ[j][2] = 0;

    for (i = 0; i < noOfPointSets; i++) {
      meanXYZ[j][0] += data[i][j][0];
      meanXYZ[j][1] += data[i][j][1];
      meanXYZ[j][2] += data[i][j][2];
    }

    meanXYZ[j][0] /= noOfPointSets;
    meanXYZ[j][1] /= noOfPointSets;
    meanXYZ[j][2] /= noOfPointSets;
  }

  // A set of 3x3 spatial covariance matrices, one for each landmark
  irtkMatrix *matrices = new irtkMatrix[nPts];
  double covarElement;

  // For each landmark:
  for (j = 0; j < nPts; j++) {
    matrices[j].Initialize(3, 3);

    for (a = 0; a < 3; a++) {
      for (b = a; b < 3; b++) {
        covarElement = 0.0;
        // Loop over sets
        for (i = 0; i < noOfPointSets; i++) {
          covarElement += (data[i][j][a] - meanXYZ[j][a])
              * (data[i][j][b] - meanXYZ[j][b]);
        }
        covarElement /= (noOfPointSets - 1.0);
        matrices[j](a, b) = covarElement;
        matrices[j](b, a) = covarElement;
      }
    }

  }

  irtkMatrix *e1 = new irtkMatrix[nPts];
  irtkMatrix *e2 = new irtkMatrix[nPts];
  irtkVector *eVals = new irtkVector[nPts];

  for (j = 0; j < nPts; j++) {
    matrices[j].Eigenvalues(e1[j], eVals[j], e2[j]);

    for (a = 0; a < 3; a++) {
      // square root of variance
      eVals[j](a) = sqrt(eVals[j](a));

      // Scale the eigenvector.
      // The scale is interpreted as a diameter by the tensor glyph.
      for (b = 0; b < 3; b++) {
        e1[j](b, a) *= (2.0 * radius * eVals[j](a));

      }
    }

  }



  // Locations for the ellipsoids, each will be centred
  // at the mean location of the corresponding landmark

  vtkPoints *outPoints = vtkPoints::New();

  outPoints->SetNumberOfPoints(nPts);
  for (j = 0; j < nPts; j++) {
    p[0] = meanXYZ[j][0];
    p[1] = meanXYZ[j][1];
    p[2] = meanXYZ[j][2];
    outPoints->InsertPoint(j, p);
  }

  // Storage for the tensor data that will be passed to the glyph filter.
  vtkDoubleArray *covTensors = vtkDoubleArray::New();
  covTensors->SetNumberOfTuples(nPts);
  covTensors->SetNumberOfComponents(9);

  for (j = 0; j < nPts; j++) {

    for (a = 0; a < 3; a++) {
      tensorVec[a]     = e1[j](a, 0);
      tensorVec[a + 3] = e1[j](a, 1);
      tensorVec[a + 6] = e1[j](a, 2);
    }

    covTensors->InsertTuple(j, tensorVec);
  }

  covTensors->Squeeze();
  covTensors->SetName("cov");

  vtkPolyData *pdWithTensors = vtkPolyData::New();
  pdWithTensors->SetPoints(outPoints);

  pdWithTensors->GetPointData()->SetTensors(covTensors);

  // The glyph to be used.
  vtkSphereSource *sphere = vtkSphereSource::New();
  sphere->SetThetaResolution(30);
  sphere->SetPhiResolution(30);

  // The glyph filter.
  vtkTensorGlyph *glyphs = vtkTensorGlyph::New();

  glyphs->SetSourceData(sphere->GetOutput());
  glyphs->SetInputData(pdWithTensors);

  // Don't do any calculation in the filter,
  // we have done all the calculations already.
  glyphs->SetExtractEigenvalues(0);
  glyphs->SetThreeGlyphs(0);
  glyphs->SetScaling(0);
  glyphs->Update();

  vtkPolyDataWriter *pdWriter = vtkPolyDataWriter::New();
  pdWriter->SetFileName(outputFileName);
  pdWriter->SetInputData(glyphs->GetOutput());
  pdWriter->SetFileTypeToBinary();
  pdWriter->Update();
  pdWriter->Write();

  Deallocate<double>(data);
  Deallocate<double>(meanXYZ);

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ) {
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
