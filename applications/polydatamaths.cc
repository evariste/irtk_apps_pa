/*
 * polydatamaths.cc
 *
 *  Created on: Jan 12, 2012
 *      Author: Paul Aljabar
 */

#ifdef HAS_VTK

#include <sys/stat.h>

#include <vtkIndent.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>


// Functions can only be binary or unary.
typedef double (*binaryFunc)(double, double);
typedef double (*unaryFunc)(double);


#define NUMERIC_OPERAND 1
#define ARRAY_OPERAND 2

char *input_name = NULL;
char *output_name = NULL;
char *current_scalars_name = NULL;
char *output_array_name = NULL;

//////////////////////////////////////////////////////////////
//
// Can add new operations if they are needed.
//
// Adding a new operation only affects Sections 1 to 3
//
//////////////////////////////////////////////////////////////

//
// Section 0: Helper functions.
//

bool is_vtkPolyDataFile(const char* filename)
{
  // Rough and ready check to see if a file contains polydata.
  char *c = NULL;
  char buff[1000];
  struct stat statStruct ;

  sprintf(buff, "%s", filename);
  c = strstr(buff, ".vtk\0");

  if (c == NULL){
    return false;
  }

  int i = stat( buff , &statStruct );
  if( i != 0 ) {
    return false;
  }


  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(filename);
  if (! reader->IsFilePolyData()){
    return false;
  }

  return true;
}

bool is_numeric(const char *str)
{
  char *pEnd;
  double dummy;
  dummy = strtod(str, &pEnd);

  if (*pEnd != '\0')
    return false;

  return true;
}



//
// Section 1 : Enum for all operations. Update number of operations when adding a new operation.
//

typedef enum {OP_NULL,
  OP_ADD,
  OP_SUB,
  OP_MUL,
  OP_DIV,
  OP_SQRT,
  OP_LOG,
  OP_ABS,
  OP_UTHR,
  OP_LTHR,
  OP_BIN
} abcdOperationType;

#define NUMBER_OF_OPERATIONS 11

//
// Section 2 : Function definitions.
//

double add_wrap(double a, double b){
  return a + b;
}

double sub_wrap(double a, double b){
  return a - b;
}

double mul_wrap(double a, double b){
  return a * b;
}

double div_wrap(double a, double b){
  return a / b;
}

double sqrt_wrap(double a){
  return sqrt(a);
}

double log_wrap(double a){
  return log(a);
}

double fabs_wrap(double a){
  return fabs(a);
}

double uthresh(double a, double thresh){
  return a > thresh ? 0.0 : a;
}

double lthresh(double a, double thresh){
  return a < thresh ? 0.0 : a;
}

double bin(double a){
  return a > 0 ? 1.0 : 0.0;
}


//
// Section 3: Functions to manipulate operations.
//


// Helper for parsing command line flags.
bool is_operation_flag(const char *str){

  if ( (strcmp(str, "-add") == 0) )
    return true;

  if ( (strcmp(str, "-sub") == 0) )
    return true;

  if ( (strcmp(str, "-mul") == 0) )
    return true;

  if ( (strcmp(str, "-div") == 0) )
    return true;

  if ( (strcmp(str, "-log") == 0) )
    return true;

  if ( (strcmp(str, "-sqrt") == 0) )
    return true;

  if ( (strcmp(str, "-abs") == 0) )
    return true;

  if ( (strcmp(str, "-uthr") == 0) )
    return true;

  if ( (strcmp(str, "-lthr") == 0) )
    return true;

  if ( (strcmp(str, "-bin") == 0) )
    return true;

  return false;
}

// Return enum identifier of an operation given its string.
abcdOperationType get_operation(const char *str){

  if ( (strcmp(str, "-add") == 0) )
    return OP_ADD;

  if ( (strcmp(str, "-sub") == 0) )
    return OP_SUB;

  if ( (strcmp(str, "-mul") == 0) )
    return OP_MUL;

  if ( (strcmp(str, "-div") == 0) )
    return OP_DIV;

  if ( (strcmp(str, "-log") == 0) )
    return OP_LOG;

  if ( (strcmp(str, "-sqrt") == 0) )
    return OP_SQRT;

  if ( (strcmp(str, "-abs") == 0) )
    return OP_ABS;

  if ( (strcmp(str, "-uthr") == 0) )
    return OP_UTHR;

  if ( (strcmp(str, "-lthr") == 0) )
    return OP_LTHR;

  if ( (strcmp(str, "-bin") == 0) )
    return OP_BIN;

  // If we got to here, there is a problem.
  return OP_NULL;
}


// Print an operation's string description given its enum position.
void print_operation(ostream& os, abcdOperationType op)
{
  switch (op){
  case OP_ADD:
    os << "add   ";
    break;
  case OP_SUB:
    os << "sub   ";
    break;
  case OP_MUL:
    os << "mul   ";
    break;
  case OP_DIV:
    os << "div   ";
    break;
  case OP_LOG:
    os << "log   ";
    break;
  case OP_SQRT:
    os << "sqrt  ";
    break;
  case OP_ABS:
    os << "abs   ";
    break;
  case OP_UTHR:
    os << "uthr  ";
    break;
  case OP_LTHR:
    os << "lthr  ";
    break;
  case OP_BIN:
    os << "bin   ";
    break;
  default:
    os << "unknown";
    break;
  }
}

// Give some help text for an operation.  First few are not
// really needed so just print a blank string.
void print_operation_help(ostream& os, abcdOperationType op)
{
  switch (op){
  case OP_ADD:
  case OP_SUB:
  case OP_MUL:
  case OP_DIV:
  case OP_LOG:
  case OP_SQRT:
    os << "  ";
    break;
  case OP_ABS:
    os << "Absolute value   ";
    break;
  case OP_UTHR:
    os << "Upper threshold, set anything greater than next operand to zero.";
    break;
  case OP_LTHR:
    os << "Lower threshold, set anything lower than next operand to zero.";
    break;
  case OP_BIN:
    os << "Binarise, all positive values are set to one. Otherwise zero.   ";
    break;
  default:
    os << "unknown";
    break;
  }
}

bool is_binary_operation(abcdOperationType op){
  switch (op){
  case OP_ADD:
  case OP_SUB:
  case OP_MUL:
  case OP_DIV:
  case OP_UTHR:
  case OP_LTHR:
    return true;
    break;
  default:
    return false;
  }
}

bool is_unary_operation(abcdOperationType op){
  switch(op){
  case OP_LOG:
  case OP_SQRT:
  case OP_ABS:
  case OP_BIN:
    return true;
    break;
  default:
    return false;
  }
}

// Return a function pointer given the enum position (Binary).
binaryFunc get_binary_operation(abcdOperationType op)
{
  binaryFunc f;

  switch (op){
  case OP_ADD:
    f = &add_wrap;
    break;
  case OP_SUB:
    f = &sub_wrap;
    break;
  case OP_MUL:
    f = &mul_wrap;
    break;
  case OP_DIV:
    f = div_wrap;
    break;
  case OP_UTHR:
    f = uthresh;
    break;
  case OP_LTHR:
    f = lthresh;
    break;
  default:
    cerr << "Invalid binary operation" << endl;
    exit(1);
  }
  return f;
}

// Return a function pointer given the enum position (Unary).
unaryFunc get_unary_operation(abcdOperationType op){
  unaryFunc f;

  switch (op){
  case OP_LOG:
    f = &log_wrap;
    break;
  case OP_SQRT:
    f = &sqrt_wrap;
    break;
  case OP_ABS:
    f = &fabs_wrap;
    break;
  case OP_BIN:
    f = bin;
    break;
  default:
    cerr << "Invalid unary operation" << endl;
    exit(1);
  }

  return f;
}


//
// Section 4: Functions to actually execute operations.
//


// Binary operation on two arrays: Result in arr1.
// Pre: both arrays are single component arrays.
void apply_operation(vtkFloatArray* arr1, vtkFloatArray *arr2, abcdOperationType op, bool &nanVals, double nanReplacement = 0)
{
  int i, noOfPoints;
  binaryFunc binFunc = NULL;
  double val;

  binFunc = get_binary_operation(op);
  nanVals = false;

  noOfPoints = arr1->GetNumberOfTuples();
  for (i = 0; i < noOfPoints; ++i){
    val = binFunc(arr1->GetTuple1(i), arr2->GetTuple1(i));

    if (isnan(val)){
      nanVals	= true;
      val = nanReplacement;
    }
    arr1->SetTuple1(i, val);

  }
}

// Binary operation on an array and a single scalar: Result in arr1.
void apply_operation(vtkFloatArray* arr, float fVal, abcdOperationType op, bool &nanVals, double nanReplacement = 0)
{
  int i, noOfPoints;
  binaryFunc binFunc;
  double val;

  binFunc = get_binary_operation(op);
  nanVals = false;

  noOfPoints = arr->GetNumberOfTuples();
  for (i = 0; i < noOfPoints; ++i){
    val = binFunc(arr->GetTuple1(i), fVal);

    if (isnan(val)){
      nanVals	= true;
      val = nanReplacement;
    }
    arr->SetTuple1(i, val);
  }
}

// Unary operation: Result in arr1.
void apply_operation(vtkFloatArray* arr, abcdOperationType op, bool &nanVals, double nanReplacement = 0)
{
  int i, noOfPoints;
  unaryFunc uniFunc;
  double val;

  uniFunc = get_unary_operation(op);
  nanVals = false;

  noOfPoints = arr->GetNumberOfTuples();
  for (i = 0; i < noOfPoints; i++){
    val = uniFunc(arr->GetTuple1(i));
    if (isnan(val)){
      nanVals	= true;
      val = nanReplacement;
    }
    arr->SetTuple1(i, val);

  }
}

//
// Section 5: And the rest ...
//

void usage()
{
  cerr << "Usage:  polydatamaths [polydataIN] [op1] [arg1] <[op2] [arg2]> ... <option(s)> ...  [polydataOUT]" << endl;
  cerr << "" << endl;
  cerr << "Apply a set of arithmetic operations to scalars associated with a polydata." << endl;
  cerr << "" << endl;
  cerr << "First argument must be a polydata file." << endl;
  cerr << "Last argument is the name of a polydata file to write the result to." << endl;
  cerr << "" << endl;
  cerr << "Other flags:" << endl;
  cerr << "    -array_name [s]        Use the array with name 's'. Applies to the last polydata set read in." << endl;
  cerr << "    -output_array_name [s] Give the final output array the name 's'" << endl;
  cerr << "    -repForNan [value]     Replacement value to use if calculations result in NAN." << endl;
  cerr << "    -q                     Quiet" << endl;
  cerr << "" << endl;
  cerr << "Possible Operations:" << endl;

  // Start from 1 so that we skip OP_NULL.
  for (int i = 1; i < NUMBER_OF_OPERATIONS; i++){
    cerr << "    -";
    print_operation(std::cerr, (abcdOperationType)i);
    print_operation_help(std::cerr, (abcdOperationType)i);
    cerr << endl;
  }
  cerr << "" << endl;

  exit(1);
}


int main(int argc, char **argv ){

  if (argc < 3 ){
    usage();
  }

  int i, arrayCount, arrayIndex, scalarCount, scalarIndex, operandIndex, opCount, opIndex;
  double nanReplacement(0.0f);
  int iVal;
  bool ok;
  int noOfPoints;
  bool nanVals;
  bool verbose = true;

  ///////////////////////////////////////////////
  //
  // Initial parsing of command line arguments.
  //


  // How many and what type of arguments do we have?
  // A pre-reading check of all arguments except last (output)
  // which must be a polydata file (but may not yet exist).

  arrayCount  = 0;
  scalarCount = 0;
  opCount     = 0;

  for (i = 1; i < argc - 1; i++){

    if (is_vtkPolyDataFile(argv[i])){
      arrayCount++;
    }

    if (is_numeric(argv[i])){
      if (strcmp(argv[i-1], "-repForNan") == 0){
        // Ignore this number, it is replacement value for nans
      } else {
        scalarCount++;
      }
    }

    if (strcmp(argv[i], "-q") == 0){
      // Deal with this argument here rather than in main part because it affects output to terminal.
      cout << "verbose to false" << endl;
      verbose = false;
    }

    if (is_operation_flag(argv[i]))
      opCount++;
  }


  if (verbose){
    cout << "Arrays  : " << arrayCount << endl;
    cout << "Scalars : " << scalarCount << endl;
    cout << "Ops     : " << opCount << endl;
  }


  vtkPolyData **pds      = new vtkPolyData*[arrayCount];

  vtkFloatArray **arrays = new vtkFloatArray* [arrayCount];
  double *scalars        = new double[scalarCount];
  int *operandTypes      = new int[scalarCount + arrayCount];

  abcdOperationType *ops = new abcdOperationType[opCount];

  for (i = 0; i < arrayCount; i++){
    pds[i]    = NULL;
    arrays[i] = vtkFloatArray::New();
  }

  ///////////////////////////////////////////////
  //
  // Deal with first argument which must be a polydata file.
  //

  input_name = argv[1];
  argv++;
  argc--;

  if (!is_vtkPolyDataFile(input_name)){
    cerr << "\n\tpolydatamaths: First argument must be a polydata file" << endl << endl;
    usage();
  }

  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(input_name);
  pds[0] = vtkPolyData::New();
  pds[0] = reader->GetOutput();
  pds[0]->Update();

  noOfPoints = pds[0]->GetNumberOfPoints();

  if (pds[0]->GetPointData()->GetNumberOfArrays() < 1){
    cerr << endl;
    cerr << "  No arrays in first data set: " << input_name << endl;
    cerr << "  Creating default set with name 'scalars'" << endl;
    arrays[0] = vtkFloatArray::New();
    arrays[0]->SetNumberOfComponents(1);
    arrays[0]->SetNumberOfTuples(noOfPoints);
    arrays[0]->SetName("scalars");

  } else {
    // Read in first set of scalars.
    arrays[0]->DeepCopy(pds[0]->GetPointData()->GetArray(0));
    arrays[0]->SetName(pds[0]->GetPointData()->GetArray(0)->GetName());
  }

  operandTypes[0] = ARRAY_OPERAND;


  ///////////////////////////////////////////////
  //
  // Main parsing of the command line arguments:
  //

  arrayIndex   =  0;
  scalarIndex  = -1;
  operandIndex =  0;
  opIndex      = -1;


  while (argc > 2){

    ok = false;

    // Do we have a polydata file?
    if ((ok == false) && is_vtkPolyDataFile(argv[1]) ){
      arrayIndex++;
      operandIndex++;

      input_name = argv[1];
      argc--;
      argv++;

      if (verbose){
        cout << "Reading polydata : " << input_name << endl;
      }

      reader->SetFileName(input_name);
      pds[arrayIndex] = vtkPolyData::New();
      pds[arrayIndex] = reader->GetOutput();
      pds[arrayIndex]->Update();

      if (pds[arrayIndex]->GetNumberOfPoints() != noOfPoints){
        cerr << "\n\tAll polydata must have the same number of points each." << endl << endl;
        exit(1);
      }
      if (pds[arrayIndex]->GetPointData()->GetNumberOfArrays() < 1){
        cerr << "\n\t No arrays in polydata " << input_name << endl << endl;
        exit(1);
      }

      vtkFloatArray *temp = static_cast<vtkFloatArray*>
      (pds[arrayIndex]->GetPointData()->GetArray(0));

      if (temp == NULL){
        cerr << "\n\tScalars unavailable";
        cerr << " in polydata object " << input_name << endl << endl;
        usage();
      }

      arrays[arrayIndex]->DeepCopy(temp);
      arrays[arrayIndex]->SetName(temp->GetName());
      operandTypes[operandIndex] = ARRAY_OPERAND;
      ok = true;
    }


    // Do we have a number?
    if ((ok == false) && is_numeric(argv[1])){
      scalarIndex++;
      operandIndex++;

      scalars[scalarIndex] = atof(argv[1]);
      argc--;
      argv++;
      operandTypes[operandIndex] = NUMERIC_OPERAND;

      ok = true;
    }

    // Do we have an operation?
    if ((ok == false) && is_operation_flag(argv[1])){
      opIndex++;
      ops[opIndex] = get_operation(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    // Do we have a name for an array?
    if ((ok == false) && (strcmp(argv[1], "-array_name") == 0)){
      // Named set of scalars - applies to last read polydata set.
      argc--;
      argv++;

      current_scalars_name = argv[1];
      argc--;
      argv++;

      vtkFloatArray *temp = static_cast<vtkFloatArray*>
      (pds[arrayIndex]->GetPointData()->GetArray(current_scalars_name, iVal));

      if (temp == NULL || iVal == -1){
        cerr << "\n\tScalars unavailable with name ";
        cerr << argv[1] << " in polydata object " << input_name << endl << endl;
        exit(1);
      }

      arrays[arrayIndex]->DeepCopy( temp );
      arrays[arrayIndex]->SetName(temp->GetName());

      ok = true;
    }


    // Name for output scalar array.
    if ((ok == false) && (strcmp(argv[1], "-output_array_name") == 0)){
      argc--;
      argv++;

      output_array_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }

    // Numeric value to replace an NAN values that occur during calculation.
    if ((ok == false) && (strcmp(argv[1], "-repForNan") == 0)){
      argc--;
      argv++;

      nanReplacement = atof(argv[1]);
      cerr << "nan replacement value = " << nanReplacement << endl;
      argc--;
      argv++;
      ok = true;
    }

    // Suppress output
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      // Already dealt with above.
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Last argument is name of output file.
  output_name = argv[1];
  argv++;
  argc--;


  ///////////////////////////////////////////////

  arrayIndex = -1;
  scalarIndex = -1;

  cout << endl << "Operands: " << endl;

  for (i = 0; i < scalarCount+arrayCount; i++){
    switch (operandTypes[i]){
    case ARRAY_OPERAND:
      arrayIndex++;
      if (verbose){
        cout << "Array   : " << arrays[arrayIndex]->GetName() << endl;
      }
      break;
    case NUMERIC_OPERAND:
      scalarIndex++;
      if (verbose){
        cout << "Numeric : " << scalars[scalarIndex] << endl;
      }
      break;
    default:
      cerr << "Error   : invalid operand type" << endl;
      exit(1);
    }
  }

  cout << endl;

  // TODO: check all arrays are only single component arrays.

  ///////////////////////////////////////////////


  if (verbose)
    cout << "Applying operations: " << endl;

  vtkFloatArray *arrayOut = vtkFloatArray::New();
  arrayOut->DeepCopy(arrays[0]);

  arrayIndex   =  0;
  scalarIndex  = -1;
  operandIndex =  0;
  opIndex      = -1;
  abcdOperationType currOp;

  // Main bit
  while (opIndex < opCount - 1){
    // Get next operation
    opIndex++;
    currOp = ops[opIndex];

    if (verbose){
      cout << "Op.     : ";
      print_operation(cout, currOp);
      cout << endl;
    }

    nanVals = false;

    if (is_binary_operation(currOp)){
      // Get next operand
      operandIndex++;
      if (operandIndex >= scalarCount+arrayCount){
        cerr << "Operand index exceeds limit: run out of operands" << endl;
        exit(1);
      }

      // Apply operation and store result in out array
      switch (operandTypes[operandIndex]){
      case ARRAY_OPERAND:
        arrayIndex++;
        if (arrayIndex >= arrayCount){
          cerr << "Array index exceeds limit: run out of array operands" << endl;
          exit(1);
        }
        apply_operation(arrayOut, arrays[arrayIndex], currOp, nanVals, nanReplacement);
        break;
      case NUMERIC_OPERAND:
        scalarIndex++;
        if (scalarIndex >= scalarCount){
          cerr << "Scalar value index exceeds limit: run out of scalar operands" << endl;
          exit(1);
        }
        apply_operation(arrayOut, scalars[scalarIndex], currOp, nanVals, nanReplacement);
        break;
      default:
        cerr << "Invalid binary operand type." << endl;
        exit(1);
      }

    } else if (is_unary_operation(currOp)){
      // Apply unary operation and modify arrayOut in place.
      apply_operation(arrayOut, currOp, nanVals, nanReplacement);

    } else {
      cerr << "Operations may only be binary or unary." << endl;
      exit(1);
    }

    if (nanVals && verbose){
      cerr << "Nan values applying current operation " << endl;
      cerr << "replaced with " << nanReplacement << endl;
    }

  }


  ///////////////////////////////////////////////


  if (verbose){
    cout << endl;
    cout << "Output  : " << output_name << endl;
  }

  // Include the resulting array and write output

  // Clear existing arrays out.
  while(pds[0]->GetPointData()->GetNumberOfArrays() > 0){
    vtkFloatArray *currArray;
    currArray = (vtkFloatArray*) pds[0]->GetPointData()->GetArray(0);
    pds[0]->GetPointData()->RemoveArray(currArray->GetName());
  }

  // Name and append the output array.
  if (output_array_name == NULL){
    arrayOut->SetName(arrays[0]->GetName());
  } else {
    arrayOut->SetName(output_array_name);
  }

  pds[0]->GetPointData()->AddArray(arrayOut);
  pds[0]->Modified();
  pds[0]->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(output_name);
  writer->SetInput(pds[0]);
  writer->Write();

  ///////////////////////////////////////////////

  delete [] pds;
  delete [] arrays;

  exit(0);

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
