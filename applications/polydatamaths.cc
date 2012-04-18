/*
 * polydatamaths.cc
 *
 *  Created on: Jan 12, 2012
 *      Author: paul
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

void usage(){
  cerr << "polydatamaths [input] op operand op operand ... [output]"<<endl;
  cerr << "Options:" << endl;
  exit(1);
}

typedef double (*binaryFunc)(double, double);
typedef double (*unaryFunc)(double);


#define NUMERIC_OPERAND 1
#define ARRAY_OPERAND 2

typedef enum {OP_NULL,
							OP_ADD,
							OP_SUB,
							OP_MUL,
							OP_DIV,
							OP_SQRT,
							OP_LOG,
							OP_ABS
						 } abcdOperationType;


char *input_name = NULL;
char *output_name = NULL;
char *current_scalars_name = NULL;
char *output_array_name = NULL;

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

bool is_operation(const char *str){

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

	return false;
}

bool is_binary_operation(abcdOperationType op){
	switch (op){
	case OP_ADD:
	case OP_SUB:
	case OP_MUL:
	case OP_DIV:
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
		return true;
		break;
	default:
		return false;
	}
}

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

	// If we got to here, there is a problem.
	return OP_NULL;
}

void print_operation(ostream& os, abcdOperationType op)
{
	switch (op){
	case OP_ADD:
		os << "add";
		break;
	case OP_SUB:
		os << "sub";
		break;
	case OP_MUL:
		os << "mul";
		break;
	case OP_DIV:
		os << "div";
		break;
	case OP_LOG:
		os << "log";
		break;
	case OP_SQRT:
		os << "sqrt";
		break;
	case OP_ABS:
		os << "abs";
		break;
	default:
		os << "unknown";
		break;
	}
}

binaryFunc get_binary_function(abcdOperationType op)
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
	default:
		cerr << "Invalid binary operation" << endl;
		exit(1);
	}
	return f;
}

unaryFunc get_unary_function(abcdOperationType op){
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
	default:
		cerr << "Invalid unary operation" << endl;
		exit(1);
	}

	return f;
}

// Binary operation on two arrays: Result in arr1.
// Pre: both arrays are single component arrays.
void apply_operation(vtkFloatArray* arr1, vtkFloatArray *arr2, abcdOperationType op, bool &nanVals, double nanReplacement = 0)
{
	int i, noOfPoints;
	binaryFunc binFunc = NULL;
	double val;

	binFunc = get_binary_function(op);
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

	binFunc = get_binary_function(op);
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

	uniFunc = get_unary_function(op);
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

  ///////////////////////////////////////////////

  // How many and what type of arguments do we have?
  // A pre-reading check of all arguments except last (output)
  // which can be a yet to exist polydata file.

  arrayCount  = 0;
  scalarCount = 0;
  opCount     = 0;

  for (i = 1; i < argc - 1; i++){

  	if (is_vtkPolyDataFile(argv[i]))
  		arrayCount++;

  	if (is_numeric(argv[i])){
  		if (strcmp(argv[i-1], "-repForNan") == 0){
  			// Ignore this number, it is replacement value for nans
  		} else {
  			scalarCount++;
  		}
  	}

  	if (is_operation(argv[i]))
  		opCount++;
  }

  cout << "Arrays:  " << arrayCount << endl;
  cout << "Scalars: " << scalarCount << endl;
  cout << "Ops    : " << opCount << endl;

  ///////////////////////////////////////////////

  vtkPolyData **pds      = new vtkPolyData*[arrayCount];

  vtkFloatArray **arrays = new vtkFloatArray* [arrayCount];
  double *scalars        = new double[scalarCount];
  int *operandTypes      = new int[scalarCount + arrayCount];

  abcdOperationType *ops = new abcdOperationType[opCount];

  for (i = 0; i < arrayCount; i++){
  	pds[i]    = NULL;
    arrays[i] = vtkFloatArray::New();
  }

  // First input has to be a polydata file
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
  if (pds[0]->GetPointData()->GetNumberOfArrays() < 1){
  	cerr << "\n\t No arrays in polydata " << input_name << endl << endl;
  	exit(1);
  }


  noOfPoints = pds[0]->GetNumberOfPoints();

  // Read in first set of scalars.
  arrays[0]->DeepCopy(pds[0]->GetPointData()->GetArray(0));
  arrays[0]->SetName(pds[0]->GetPointData()->GetArray(0)->GetName());

  operandTypes[0] = ARRAY_OPERAND;

  arrayIndex   =  0;
  scalarIndex  = -1;
  operandIndex =  0;
  opIndex      = -1;

  ///////////////////////////////////////////////

  while (argc > 2){
  	cout << argv[1] << endl;

    ok = false;

    // Do we have polydata?
    if ((ok == false) && is_vtkPolyDataFile(argv[1]) ){
      arrayIndex++;
      operandIndex++;

      input_name = argv[1];
      argc--;
      argv++;

      cout << "Reading polydata : " << input_name << endl;
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

    // Do we have a name for an array?
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)){
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

    // Do we have an operation?
    if ((ok == false) && is_operation(argv[1])){

      opIndex++;
      ops[opIndex] = get_operation(argv[1]);
      argc--;
      argv++;
    	ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-output_array_name") == 0)){
      argc--;
      argv++;

    	output_array_name = argv[1];
      argc--;
      argv++;
    	ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-repForNan") == 0)){
      argc--;
      argv++;

      nanReplacement = atof(argv[1]);
      cerr << "nan replacement value = " << nanReplacement << endl;
      argc--;
      argv++;
    	ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  ///////////////////////////////////////////////

  arrayIndex = -1;
  scalarIndex = -1;

  for (i = 0; i < scalarCount+arrayCount; i++){
  	switch (operandTypes[i]){
  	case ARRAY_OPERAND:
  		arrayIndex++;
  		cout << "Array : " << endl;
  		cout << arrays[arrayIndex]->GetName() << endl;
  		cout << "....." << endl;
  		break;
  	case NUMERIC_OPERAND:
  		scalarIndex++;
  		cout << "Numeric : " << endl;
  		cout << scalars[scalarIndex] << endl;
  		cout << "....." << endl;
  		break;
  	default:
  		cerr << "Error : invalid operand type" << endl;
  		exit(1);
  	}
  }

  // TODO: check all arrays are only single component arrays.

  ///////////////////////////////////////////////


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

  	cout << "current operation ";
  	print_operation(cout, currOp);
  	cout << endl;

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

  	if (nanVals){
  		cerr << "Nan values applying current operation " << endl;
  		cerr << "replaced with " << nanReplacement << endl;
  	}

  }

  ///////////////////////////////////////////////

  output_name = argv[1];
  argv++;
  argc--;

  cout << "Output: " << output_name << endl;

  // Add the resulting array and write output
	while(pds[0]->GetPointData()->GetNumberOfArrays() > 0){
		vtkFloatArray *currArray;
		currArray = (vtkFloatArray*) pds[0]->GetPointData()->GetArray(0);
		pds[0]->GetPointData()->RemoveArray(currArray->GetName());
	}

	if (output_array_name == NULL){
		arrayOut->SetName("result");
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
