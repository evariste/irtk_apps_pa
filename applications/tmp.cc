#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVRMLImporter.h>
#include <vtkDataSet.h>
#include <vtkActorCollection.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkLight.h>
 
int main ( int argc, char *argv[])
{
  if(argc != 2)
    {
    std::cout << "Required arguments: Filename" << std::endl;
    return EXIT_FAILURE;
    }
 
  std::string filename = argv[1];
  std::cout << "Reading " << filename << std::endl;
 
 vtkSmartPointer<vtkRenderer> ren1= vtkSmartPointer<vtkRenderer>::New();
 vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
 renderWindow->AddRenderer(ren1);

 vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  // VRML Import
  vtkSmartPointer<vtkVRMLImporter> importer = vtkSmartPointer<vtkVRMLImporter>::New();
  importer->SetFileName ( filename.c_str() );
  importer->Read();
  importer->SetRenderWindow(renderWindow);
  importer->Update();
 
ren1=importer->GetRenderer();
vtkRendererCollection *renCollection=vtkRendererCollection::New();
renCollection = renderWindow->GetRenderers();
renCollection->InitTraversal();
vtkRenderer *ren= vtkRenderer::New();
ren = renCollection->GetNextItem();

ren->ResetCamera();
vtkLight *light1=vtkLight::New();
light1->SetLightTypeToCameraLight ();
ren->AddLight(light1);

vtkActorCollection *actorcol=vtkActorCollection::New();
actorcol=ren->GetActors();

vtkActor *actor=vtkActor::New();
actor=actorcol->GetLastActor();

vtkMapper *map=actor->GetMapper();
vtkDataSet *PolyData=map->GetInput();

ren->AddActor(actor);
renderWindow->Render();
renderWindowInteractor->Start();
 
  return 0;
}




// #include <irtkImage.h>
// #include <irtkFileToImage.h>
// #include <nifti1_io.h>

// char *input_name = NULL;


// void nifti_mat44_to_orientation2( mat44 R , int *icod, int *jcod, int *kcod );


// void usage()
// {
//   cerr << "USAGE" << endl;
//   exit(1);
// }

// int main(int argc, char **argv)
// {
// 	int i,j;
//   if (argc < 2){
//     usage();
//   }

//   input_name = argv[1];
//   argc--;
//   argv++;


//   irtkFileToImage *reader = irtkFileToImage::New(input_name);
//   irtkBaseImage *img;
//   img = reader->GetOutput();

//   irtkImageAttributes attr;

//   attr = img->GetImageAttributes();

//   mat44 M;

//   irtkMatrix i2w = img->GetImageToWorldMatrix();

//   for (j=0; j<4; j++){
// 	  for (i=0; i<4; i++){
// 	    M.m[i][j] = i2w(i,j);
// 	  }
//   }

//   int icod, jcod, kcod;
//   nifti_mat44_to_orientation2(M, &icod, &jcod, &kcod);


//   cout << icod << " " << jcod << " " << kcod << endl;
//   cout << nifti_orientation_string(icod) <<  endl;
//   cout << nifti_orientation_string(jcod) <<  endl;
//   cout << nifti_orientation_string(kcod) <<  endl;



// }



// void nifti_mat44_to_orientation2( mat44 R , int *icod, int *jcod, int *kcod )
// {
//    float xi,xj,xk , yi,yj,yk , zi,zj,zk , val,detQ,detP ;
//    mat33 P , Q , M ;
//    int i,j,k=0,p,q,r , ibest,jbest,kbest,pbest,qbest,rbest ;
//    float vbest ;

//    if( icod == NULL || jcod == NULL || kcod == NULL ) return ; /* bad */

//    *icod = *jcod = *kcod = 0 ; /* error returns, if sh*t happens */

//    /* load column vectors for each (i,j,k) direction from matrix */

//    /*-- i axis --*/ /*-- j axis --*/ /*-- k axis --*/

//    xi = R.m[0][0] ; xj = R.m[0][1] ; xk = R.m[0][2] ;
//    yi = R.m[1][0] ; yj = R.m[1][1] ; yk = R.m[1][2] ;
//    zi = R.m[2][0] ; zj = R.m[2][1] ; zk = R.m[2][2] ;

//    /* normalize column vectors to get unit vectors along each ijk-axis */

//    /* normalize i axis */

//    val = sqrt( xi*xi + yi*yi + zi*zi ) ;
//    if( val == 0.0 ) return ;                 /* stupid input */
//    xi /= val ; yi /= val ; zi /= val ;

//    /* normalize j axis */

//    val = sqrt( xj*xj + yj*yj + zj*zj ) ;
//    if( val == 0.0 ) return ;                 /* stupid input */
//    xj /= val ; yj /= val ; zj /= val ;

//    /* orthogonalize j axis to i axis, if needed */

//    val = xi*xj + yi*yj + zi*zj ;    /* dot product between i and j */
//    if( fabs(val) > 1.e-4 ){
//      xj -= val*xi ; yj -= val*yi ; zj -= val*zi ;
//      val = sqrt( xj*xj + yj*yj + zj*zj ) ;  /* must renormalize */
//      if( val == 0.0 ) return ;              /* j was parallel to i? */
//      xj /= val ; yj /= val ; zj /= val ;
//    }

//    /* normalize k axis; if it is zero, make it the cross product i x j */

//    val = sqrt( xk*xk + yk*yk + zk*zk ) ;
//    if( val == 0.0 ){ xk = yi*zj-zi*yj; yk = zi*xj-zj*xi ; zk=xi*yj-yi*xj ; }
//    else            { xk /= val ; yk /= val ; zk /= val ; }

//    /* orthogonalize k to i */

//    val = xi*xk + yi*yk + zi*zk ;    /* dot product between i and k */
//    if( fabs(val) > 1.e-4 ){
//      xk -= val*xi ; yk -= val*yi ; zk -= val*zi ;
//      val = sqrt( xk*xk + yk*yk + zk*zk ) ;
//      if( val == 0.0 ) return ;      /* bad */
//      xk /= val ; yk /= val ; zk /= val ;
//    }

//    /* orthogonalize k to j */

//    val = xj*xk + yj*yk + zj*zk ;    /* dot product between j and k */
//    if( fabs(val) > 1.e-4 ){
//      xk -= val*xj ; yk -= val*yj ; zk -= val*zj ;
//      val = sqrt( xk*xk + yk*yk + zk*zk ) ;
//      if( val == 0.0 ) return ;      /* bad */
//      xk /= val ; yk /= val ; zk /= val ;
//    }

//    Q.m[0][0] = xi ; Q.m[0][1] = xj ; Q.m[0][2] = xk ;
//    Q.m[1][0] = yi ; Q.m[1][1] = yj ; Q.m[1][2] = yk ;
//    Q.m[2][0] = zi ; Q.m[2][1] = zj ; Q.m[2][2] = zk ;

//    /* at this point, Q is the rotation matrix from the (i,j,k) to (x,y,z) axes */

//    detQ = nifti_mat33_determ( Q ) ;
//    if( detQ == 0.0 ) return ; /* shouldn't happen unless user is a DUFIS */

//    /* Build and test all possible +1/-1 coordinate permutation matrices P;
//       then find the P such that the rotation matrix M=PQ is closest to the
//       identity, in the sense of M having the smallest total rotation angle. */

//    /* Despite the formidable looking 6 nested loops, there are
//       only 3*3*3*2*2*2 = 216 passes, which will run very quickly. */

//    vbest = -666.0 ; ibest=pbest=qbest=rbest=1 ; jbest=2 ; kbest=3 ;
//    for( i=1 ; i <= 3 ; i++ ){     /* i = column number to use for row #1 */
//     for( j=1 ; j <= 3 ; j++ ){    /* j = column number to use for row #2 */
//      if( i == j ) continue ;
//       for( k=1 ; k <= 3 ; k++ ){  /* k = column number to use for row #3 */
//        if( i == k || j == k ) continue ;
//        P.m[0][0] = P.m[0][1] = P.m[0][2] =
//         P.m[1][0] = P.m[1][1] = P.m[1][2] =
//          P.m[2][0] = P.m[2][1] = P.m[2][2] = 0.0 ;
//        for( p=-1 ; p <= 1 ; p+=2 ){    /* p,q,r are -1 or +1      */
//         for( q=-1 ; q <= 1 ; q+=2 ){   /* and go into rows #1,2,3 */
//          for( r=-1 ; r <= 1 ; r+=2 ){
//            P.m[0][i-1] = p ; P.m[1][j-1] = q ; P.m[2][k-1] = r ;
//            detP = nifti_mat33_determ(P) ;           /* sign of permutation */
//            if( detP * detQ <= 0.0 ) continue ;  /* doesn't match sign of Q */
//            M = nifti_mat33_mul(P,Q) ;

//            /* angle of M rotation = 2.0*acos(0.5*sqrt(1.0+trace(M)))       */
//            /* we want largest trace(M) == smallest angle == M nearest to I */

//            val = M.m[0][0] + M.m[1][1] + M.m[2][2] ; /* trace */
//            if( val > vbest ){
//              vbest = val ;
//              ibest = i ; jbest = j ; kbest = k ;
//              pbest = p ; qbest = q ; rbest = r ;
//            }
//    }}}}}}

//    P.m[0][0] = P.m[0][1] = P.m[0][2] =
//     P.m[1][0] = P.m[1][1] = P.m[1][2] =
//      P.m[2][0] = P.m[2][1] = P.m[2][2] = 0.0 ;

//    P.m[0][ibest-1] = pbest ;
//    P.m[1][jbest-1] = qbest ;
//    P.m[2][kbest-1] = rbest ;


//    cout << "====================" << endl;
//    cout << endl;

//    for (i = 0; i < 3; i++){
//      for (j  = 0; j < 3; j++){
//        cout << "   " << P.m[i][j];
//      }
//      cout << endl;
//    }

//    cout << endl;
//    cout << "====================" << endl;


//    /* At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.

//       The matrix P that corresponds is the best permutation approximation
//       to Q-inverse; that is, P (approximately) takes (x,y,z) coordinates
//       to the (i,j,k) axes.

//       For example, the first row of P (which contains pbest in column ibest)
//       determines the way the i axis points relative to the anatomical
//       (x,y,z) axes.  If ibest is 2, then the i axis is along the y axis,
//       which is direction P2A (if pbest > 0) or A2P (if pbest < 0).

//       So, using ibest and pbest, we can assign the output code for
//       the i axis.  Mutatis mutandis for the j and k axes, of course. */

//    switch( ibest*pbest ){
//      case  1: i = NIFTI_L2R ; break ;
//      case -1: i = NIFTI_R2L ; break ;
//      case  2: i = NIFTI_P2A ; break ;
//      case -2: i = NIFTI_A2P ; break ;
//      case  3: i = NIFTI_I2S ; break ;
//      case -3: i = NIFTI_S2I ; break ;
//    }

//    switch( jbest*qbest ){
//      case  1: j = NIFTI_L2R ; break ;
//      case -1: j = NIFTI_R2L ; break ;
//      case  2: j = NIFTI_P2A ; break ;
//      case -2: j = NIFTI_A2P ; break ;
//      case  3: j = NIFTI_I2S ; break ;
//      case -3: j = NIFTI_S2I ; break ;
//    }

//    switch( kbest*rbest ){
//      case  1: k = NIFTI_L2R ; break ;
//      case -1: k = NIFTI_R2L ; break ;
//      case  2: k = NIFTI_P2A ; break ;
//      case -2: k = NIFTI_A2P ; break ;
//      case  3: k = NIFTI_I2S ; break ;
//      case -3: k = NIFTI_S2I ; break ;
//    }

//    *icod = i ; *jcod = j ; *kcod = k ; return ;
// }

