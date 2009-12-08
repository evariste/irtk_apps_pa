// Trying to re-orient a ffd to two reference normalised versions
// of a target and source image.

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkBSplineFreeFormTransformation3D.h>
#include <irtkMultiLevelFreeFormTransformation.h>

// Input transformations
char *mffd_in_name = NULL;
char *aregRef2T_name = NULL;
char *aregRef2S_name = NULL;

// Output transformation
char *mffd_out_name = NULL;

void matTimesVec3by3(irtkMatrix m, double *in, double *out)
{
 // NO bounds checking!!
 int i, j;
 for (i = 0; i < 3; i++){
   out[i] = 0;
   for (j = 0; j < 3; j++){
     out[i] += m(i, j) * in[j];
   }
 }

}

void usage()
{
 cerr << "Usage: ffdreorient [mffdT2S] [aregRef2T] [aregRef2S] [mffd_out]\n" << endl;
 exit(1);
}

int main(int argc, char **argv)
{
 int xdim, ydim, zdim, i, j, k;
 double x, y, z;
 double xAffCorr, yAffCorr, zAffCorr;
 double x1, x2, y1, y2, z1, z2, dx, dy, dz;
 int xdim_tr, ydim_tr, zdim_tr;
 double x1_tr, x2_tr, y1_tr, y2_tr, z1_tr, z2_tr, dx_tr, dy_tr, dz_tr;
 double xaxis[3], yaxis[3], zaxis[3];
 double xaxis_tr[3], yaxis_tr[3], zaxis_tr[3];
 char buf[256];
 double in[3], out[3];
 double factor;


 // Check command line
 if (argc != 5){
   usage();
 }

 // Parse file names
 mffd_in_name  = argv[1];
 argc--;
 argv++;
 aregRef2T_name = argv[1];
 argc--;
 argv++;
 aregRef2S_name = argv[1];
 argc--;
 argv++;
 mffd_out_name = argv[1];
 argc--;
 argv++;

 // Read transformations
 cout << "Reading transformations ... "; cout.flush();
 irtkMultiLevelFreeFormTransformation *mffdIn = new irtkMultiLevelFreeFormTransformation;
 mffdIn->irtkTransformation::Read(mffd_in_name);

 irtkAffineTransformation *affineRef2Tgt = new irtkAffineTransformation;
 affineRef2Tgt->irtkTransformation::Read(aregRef2T_name);

 irtkAffineTransformation *affineRef2Src = new irtkAffineTransformation;
 affineRef2Src->irtkTransformation::Read(aregRef2S_name);
 cout << "done" << endl;

 if (mffdIn->NumberOfLevels() > 1){
   cerr << "Not yet implemented for MFFDs with multiple levels." << endl;
   usage();
 }

 irtkMatrix A = affineRef2Tgt->GetMatrix();
 irtkMatrix B = mffdIn->GetMatrix();
 irtkMatrix C = affineRef2Src->GetMatrix();

 irtkMatrix invC(4,4);
 invC = C;
 invC.Invert();


 irtkMatrix M = invC * B * A;

 M.Print();

 irtkRigidTransformation *rigidRef2Tgt = new irtkRigidTransformation;
 rigidRef2Tgt->PutTranslationX(affineRef2Tgt->GetTranslationX());
 rigidRef2Tgt->PutTranslationY(affineRef2Tgt->GetTranslationY());
 rigidRef2Tgt->PutTranslationZ(affineRef2Tgt->GetTranslationZ());
 rigidRef2Tgt->PutRotationX(affineRef2Tgt->GetRotationX());
 rigidRef2Tgt->PutRotationY(affineRef2Tgt->GetRotationY());
 rigidRef2Tgt->PutRotationZ(affineRef2Tgt->GetRotationZ());
 rigidRef2Tgt->UpdateMatrix();

 irtkMatrix Arigid = rigidRef2Tgt->GetMatrix();

 irtkMatrix invA(4,4);
 invA = A;
 invA.Invert();

 irtkMatrix invArigid = Arigid;
 invArigid.Invert();





 irtkMultiLevelFreeFormTransformation *mffdOut = new irtkMultiLevelFreeFormTransformation;
 mffdOut->PutMatrix(B);

 irtkFreeFormTransformation3D *ffd =
       dynamic_cast<irtkFreeFormTransformation3D *> (mffdIn->GetLocalTransformation(0));

 xdim = ffd->GetX();
 ydim = ffd->GetY();
 zdim = ffd->GetZ();

 x1 = 0;
 y1 = 0;
 z1 = 0;
 ffd->LatticeToWorld(x1, y1, z1);
 x2 = xdim - 1;
 y2 = ydim - 1;
 z2 = zdim - 1;
 ffd->LatticeToWorld(x2, y2, z2);

 ffd->GetSpacing(dx, dy, dz);

 ffd->GetOrientation(xaxis, yaxis, zaxis);

 cout << "Dims " << xdim << ", " << ydim << ", " << zdim << endl;
 cout << "Bounds (" << x1 << ", " << y1 << ", " << z1  << ") (" << x2 << ", " << y2 << ", " << z2 << ")" << endl;
 cout << "Spacing : " << dx << ", " << dy << ", " << dz << endl;
 cout << "Axes" << endl;

 sprintf(buf, "     %2.4f %2.4f %2.4f\n", xaxis[0], xaxis[1], xaxis[2]);
 cout << buf << endl;
 sprintf(buf, "     %2.4f %2.4f %2.4f\n", yaxis[0], yaxis[1], yaxis[2]);
 cout << buf << endl;
 sprintf(buf, "     %2.4f %2.4f %2.4f\n", zaxis[0], zaxis[1], zaxis[2]);
 cout << buf << endl;

 matTimesVec3by3(invArigid, xaxis, xaxis_tr);
 matTimesVec3by3(invArigid, yaxis, yaxis_tr);
 matTimesVec3by3(invArigid, zaxis, zaxis_tr);


 in[0] = x1;
 in[1] = y1;
 in[2] = z1;
 matTimesVec3by3(invA, in, out);
 x1_tr = out[0];
 y1_tr = out[1];
 z1_tr = out[2];

 in[0] = x2;
 in[1] = y2;
 in[2] = z2;
 matTimesVec3by3(invA, in, out);
 x2_tr = out[0];
 y2_tr = out[1];
 z2_tr = out[2];

 double a1, a2, b1, b2, c1, c2;

 a1 = x1_tr * xaxis_tr[0] + y1_tr * xaxis_tr[1] + z1_tr * xaxis_tr[2];
 b1 = x1_tr * yaxis_tr[0] + y1_tr * yaxis_tr[1] + z1_tr * yaxis_tr[2];
 c1 = x1_tr * zaxis_tr[0] + y1_tr * zaxis_tr[1] + z1_tr * zaxis_tr[2];

 a2 = x2_tr * xaxis_tr[0] + y2_tr * xaxis_tr[1] + z2_tr * xaxis_tr[2];
 b2 = x2_tr * yaxis_tr[0] + y2_tr * yaxis_tr[1] + z2_tr * yaxis_tr[2];
 c2 = x2_tr * zaxis_tr[0] + y2_tr * zaxis_tr[1] + z2_tr * zaxis_tr[2];


 dx_tr = (a2 - a1) / (float(xdim) - 1);
 dy_tr = (b2 - b1) / (float(ydim) - 1);
 dz_tr = (c2 - c1) / (float(zdim) - 1);


 cout << "\n before setting \n" << endl;

 cout << "Bounds (" << x1_tr << ", " << y1_tr << ", " << z1_tr  << ") (" << x2_tr << ", " << y2_tr << ", " << z2_tr << ")" << endl;
 cout << "Spacing : " << dx_tr << ", " << dy_tr << ", " << dz_tr << endl;
 cout << "Axes" << endl;

 sprintf(buf, "     %2.4f %2.4f %2.4f\n", xaxis_tr[0], xaxis_tr[1], xaxis_tr[2]);
 cout << buf << endl;
 sprintf(buf, "     %2.4f %2.4f %2.4f\n", yaxis_tr[0], yaxis_tr[1], yaxis_tr[2]);
 cout << buf << endl;
 sprintf(buf, "     %2.4f %2.4f %2.4f\n", zaxis_tr[0], zaxis_tr[1], zaxis_tr[2]);
 cout << buf << endl;


 irtkBSplineFreeFormTransformation3D *ffdOut =
     new irtkBSplineFreeFormTransformation3D(x1_tr, y1_tr, z1_tr, x2_tr, y2_tr, z2_tr, dx_tr, dy_tr, dz_tr, xaxis_tr, yaxis_tr, zaxis_tr);



 xdim_tr = ffdOut->GetX();
 ydim_tr = ffdOut->GetY();
 zdim_tr = ffdOut->GetZ();

 x1_tr = 0;
 y1_tr = 0;
 z1_tr = 0;
 ffdOut->LatticeToWorld(x1_tr, y1_tr, z1_tr);
 x2_tr = xdim_tr - 1;
 y2_tr = ydim_tr - 1;
 z2_tr = zdim_tr - 1;
 ffdOut->LatticeToWorld(x2_tr, y2_tr, z2_tr);

 ffdOut->GetSpacing(dx_tr, dy_tr, dz_tr);

 ffdOut->GetOrientation(xaxis_tr, yaxis_tr, zaxis_tr);

 cout << "\n after setting \n" << endl;

 cout << "Dims " << xdim_tr << ", " << ydim_tr << ", " << zdim_tr << endl;
 cout << "Bounds (" << x1_tr << ", " << y1_tr << ", " << z1_tr  << ") (" << x2_tr << ", " << y2_tr << ", " << z2_tr << ")" << endl;
 cout << "Spacing : " << dx_tr << ", " << dy_tr << ", " << dz_tr << endl;
 cout << "Axes" << endl;

 sprintf(buf, "     %2.4f %2.4f %2.4f\n", xaxis_tr[0], xaxis_tr[1], xaxis_tr[2]);
 cout << buf << endl;
 sprintf(buf, "     %2.4f %2.4f %2.4f\n", yaxis_tr[0], yaxis_tr[1], yaxis_tr[2]);
 cout << buf << endl;
 sprintf(buf, "     %2.4f %2.4f %2.4f\n", zaxis_tr[0], zaxis_tr[1], zaxis_tr[2]);
 cout << buf << endl;



 int noOfDofs, noOfCPs;

 noOfDofs = ffdOut->NumberOfDOFs();

 noOfCPs = noOfDofs / 3;

 if (noOfDofs != ffd->NumberOfDOFs()){
   cerr << "Error : no of dofs mismatch" << endl;
   exit(1);
 }

 for (i = 0; i < noOfCPs; i++){
   in[0] = ffd->Get(i);
   in[1] = ffd->Get(i + noOfCPs);
   in[2] = ffd->Get(i + 2 * noOfCPs);
   matTimesVec3by3(invC, in, out);
   ffdOut->Put(i, out[0]);
   ffdOut->Put(i + noOfCPs, out[1]);
   ffdOut->Put(i + 2 * noOfCPs, out[2]);

   ffdOut->PutStatus(i, ffd->GetStatus(i));
   ffdOut->PutStatus(i+ noOfCPs, ffd->GetStatus(i+ noOfCPs));
   ffdOut->PutStatus(i + 2 * noOfCPs, ffd->GetStatus(i + 2 * noOfCPs));
 }


 mffdOut->PutMatrix(M);

 mffdOut->PushLocalTransformation(ffdOut);
 mffdOut->irtkTransformation::Write(mffd_out_name);



 exit(0);


 // `de-affine' the original ffd.
 xdim = mffdIn->GetLocalTransformation(0)->GetX();
 ydim = mffdIn->GetLocalTransformation(0)->GetY();
 zdim = mffdIn->GetLocalTransformation(0)->GetZ();

 // Correct the control point values in each of the input FFDs with
 // respect to their corresponding affine components.


 // Correcting the cp values requires the affine part to be S->T.
 irtkMatrix globalIn = mffdIn->irtkAffineTransformation::GetMatrix();
 globalIn.Invert();


 for (k = 0; k < zdim; ++k){
   for (j = 0; j < ydim; ++j){
     for (i = 0; i < xdim; ++i){
       ffd->Get(i, j, k, x, y, z);

       xAffCorr = globalIn(0, 0) * x + globalIn(0, 1) * y + globalIn(0, 2) * z;
       yAffCorr = globalIn(1, 0) * x + globalIn(1, 1) * y + globalIn(1, 2) * z;
       zAffCorr = globalIn(2, 0) * x + globalIn(2, 1) * y + globalIn(2, 2) * z;

       ffd->Put(i, j, k, xAffCorr, yAffCorr, zAffCorr);
     }
   }
 }

 // apply the new affine matrix to the ffd.
 irtkMatrix globalOut = affineRef2Tgt->GetMatrix();
 mffdIn->PutMatrix(globalOut);

 ffd = dynamic_cast<irtkFreeFormTransformation3D *> (mffdIn->GetLocalTransformation(0));

 for (k = 0; k < zdim; ++k){
   for (j = 0; j < ydim; ++j){
     for (i = 0; i < xdim; ++i){
       ffd->Get(i, j, k, x, y, z);

       xAffCorr = globalOut(0, 0) * x + globalOut(0, 1) * y + globalOut(0, 2) * z;
       yAffCorr = globalOut(1, 0) * x + globalOut(1, 1) * y + globalOut(1, 2) * z;
       zAffCorr = globalOut(2, 0) * x + globalOut(2, 1) * y + globalOut(2, 2) * z;

       ffd->Put(i, j, k, xAffCorr, yAffCorr, zAffCorr);
     }
   }
 }


 // Write transformation
 mffdIn->irtkTransformation::Write(mffd_out_name);
}


