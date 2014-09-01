#ifdef HAS_VTK

#ifndef _VTKPOLYDATASMOOTHCUSTOM_H

#define _VTKPOLYDATASMOOTHCUSTOM_H


#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkCurvatures.h>

#include <irtkImage.h>

//#include <nr.h>
#include <gsl/gsl_vector.h> /*For Vectors*/
#include <gsl/gsl_sort_vector.h>

class vtkPolyDataSmoothCustom : public irtkObject
{
  vtkPolyData *_input;

  vtkPolyData *_output;

  vtkFloatArray *_normals;

  vtkFloatArray *_distances;

  double _RelaxationFactor;

  double _CofG[3];

  int _NoOfIterations;

  bool _TrackingOn;

  float _SmoothnessThreshold;

  double HSquareRobustMean();

  double SurfaceArea();

  void UpdateCentreOfGravity();

  void ShiftAndScaleMesh(double *shift, double factor);

  double GetRadius();


public:
  vtkPolyDataSmoothCustom();

  virtual ~vtkPolyDataSmoothCustom();

  void Initialize(vtkPolyData *polydata);

  void Finalize();

  void SetInput(vtkPolyData *);

  void Run();

  vtkPolyData *GetOutput();

  virtual SetMacro(NoOfIterations, int);
  virtual GetMacro(NoOfIterations, int);
  virtual SetMacro(TrackingOn, int);
  virtual GetMacro(TrackingOn, int);
  virtual SetMacro(SmoothnessThreshold, float);
  virtual GetMacro(SmoothnessThreshold, float);
  virtual SetMacro(RelaxationFactor, double);
  virtual GetMacro(RelaxationFactor, double);
};

#endif // _VTKPOLYDATASMOOTHCUSTOM_H

#endif
