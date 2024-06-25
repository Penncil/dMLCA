#include <R.h>
#include <stdlib.h>

// Test
// void hello(int *n)
// {
//  int i;
//  for(i=0; i < *n; i++) {
//    Rprintf("Hello, world!\n");
//  }
// }

//Need to write four functions:
//1. pdf/likelihood of y
//2. w (compare to see speed)
//3. updating p
//4. First/Second order derivatives for alpha and beta.

#define MAX_CLASSES 500  // Maximum number of latent classes

// Function: ylik
//       A function to find the likelihood of an observation
//       given a vector of responses and a set of
//       choice probabilities for each latent class.
//
//       Args:
//            probs (double *):   A vector of choice probabilities
//                                over all items.
//            y (int *):          A vectorize matrix of choices
//            obs (int *):        Number of observations (rows of y)
//            items (int *):      Number of items (cols of y)
//            numChoices (int *): Vector given number of choices for
//                                each item.
//            classes (int):      Number of latent classes
//
//       Return:
//            llik (double *): likelihood vector for the observation
//
// for 1 site, return a vector of length C*n: (f(y1,p1),..,f(y1,pC),f(y2,p1),...,f(yn,pC))
void ylik(double *probs, int *y, int *obs, int *items,
          int *numChoices, int *classes, double *lik) {
  int i,j,k;
  const int citems = *items;
  const int cclasses = *classes;
  const int cobs = *obs;
  const double *firstprobs = probs;
  
  for (i=0;i<cobs;i++) {
    for (j=0;j<cclasses;j++) lik[j] = 1.0;
    probs = (double *) firstprobs;
    for (k=0;k<citems;k++) {
      for (j=0;j<cclasses;j++) {
        if (y[k]>0) lik[j] *= probs[y[k]-1];
        probs += numChoices[k]; // move pointer to next item
      }
    }
    y += citems; // move pointer to next observation
    lik += cclasses; // move pointer to next observation
  }
}

// for K site, return a vector of length C*n*K: (f(y11,p1),..,f(y11,pC),f(y12,p1),...,f(y1n,pC), ..., f(yk1,p1),...)
// probs: same as before; a vector of length C * q * 2 (binary Y)
// y: a vector of length n*q*K y_knq (y111, ..., y11q, ...,y1nq, ..., yk11,...,yknq,...)
// as.vector(t(y)), as.vector not needed in fact 
void ylik_Ksites(double *probs, int *y, int*numSite, int *obs, int *items,
                 int *numChoices, int *classes, double *lik) {
  int s,i,j,k;
  const int citems = *items; //q
  const int cclasses = *classes; //C
  const int cobs = *obs; //n
  const int csite = *numSite;
  const double *firstprobs = probs;
  
  for (s=0;s<csite;s++) {
    for (i=0;i<cobs;i++) {
      for (j=0;j<cclasses;j++) lik[j] = 1.0;
      probs = (double *) firstprobs;
      for (k=0;k<citems;k++) {
        for (j=0;j<cclasses;j++) {
          if (y[k]>0) lik[j] *= probs[y[k]-1];
          probs += numChoices[k]; // move pointer to next item
        }
      }
      y += citems; // move pointer to next observation
      lik += cclasses; // move pointer to next observation
    }
  }
}
//
//           A function to find the posterior distribution over the
//           latent classes for each observation.
//
//     Args:
//           prior (double *): Vector of prior class probs
//           probs (double *): An (N x J) matrix of choice probs
//           y (int *): An (N x J) matrix of choices
//           items (int *): Number of items
//           obs (int *): Number of observations
//           numChoices (int *): Vector given number of choices for
//                                each item.
//           classes (int):      Number of latent classes
//
//     Return:
//           post (double *): Vector of poster probs
//
void postclass_Ksites(double *prior, double *probs, int *y, int*numSite,
                      int *items, int *obs, int *numChoices,
                      int *classes, double *posterior) {
  int i,j,totalChoices,s;
  double llik[MAX_CLASSES]; 
  double denom;
  const int citems = *items;
  const int cobs = *obs;
  const int cclasses = *classes;
  int one = 1;
  const int csite = *numSite;
  
  totalChoices=0;
  for (i=0;i<citems;i++) totalChoices += numChoices[i]; //K1+...+Kj
  for (s=0;s<csite;s++){
    for (i=0;i<cobs;
    i++) {
      ylik(probs,y, (int *) &one,items,numChoices,classes,llik);
      denom = 0.0;
      for (j=0;j<cclasses;j++) denom+=prior[j]*llik[j];
      for (j=0;j<cclasses;j++) {
        posterior[j]=prior[j]*llik[j]/denom;
      }
      y+=citems; //Increment y pointer to next obs
      prior+=cclasses; // Increment proir pointer to next obs
      posterior+=cclasses;
    }
  }
}
// Function: probhat
//     A function to return updates estimates the response probs
//     within each class given the data and posterior distributions
//     Just a bunch of conditional means.
//
//  Args:
//     y (int *): Matrix of response data
//     post (double *): Matrix of posterior estimates
//     items (int *): Number of items
//     obs (int *): Number of observations
//     numChoices (int *): Vector given number of choices for
//                         each item.
//     classes (int):      Number of latent classes
//
//   Returns:
//     probhat (double *): Estimates response probs in mat form
void probhat_Ksites(int *y, int *numSite, double *post,
                    int *items, int *obs, int *numChoices,
                    int *classes, double *probhat) {
  double *denom;
  int i,j,k,cumChoices,s;
  const int citems = *items;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int csite = *numSite;
  
  int totalChoices=0;
  for (i=0;i<citems;i++) totalChoices += numChoices[i];
  for (i=0;i<(totalChoices*cclasses);i++) probhat[i]=0.0;
  
  denom = (double *) calloc((cclasses*citems),sizeof(double));
  for (i=0;i<(cclasses*citems);i++) denom[i]=0.0;
  
  for (s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      for (j=0;j<cclasses;j++) {
        cumChoices=0;
        for (k=0;k<citems;k++) {
          if (y[k]>0) {
            probhat[j*numChoices[k]+cclasses*cumChoices+y[k]-1] += post[j];
            denom[j*citems+k] += post[j];
          }
          cumChoices += numChoices[k];
        }
      }
      y+=citems;
      post+=cclasses;
    }
  }
  
  for (j=0;j<cclasses;j++) {
    cumChoices=0;
    for (k=0;k<citems;k++) {
      for (i=0;i<numChoices[k];i++) {
        probhat[j*numChoices[k]+cclasses*cumChoices+i] =
          (double) probhat[j*numChoices[k]+cclasses*cumChoices+i]/
            denom[j*citems+k];
      }
      cumChoices += numChoices[k];
    }
  }
  free(denom);
}
// Function: d2lldbeta2
//	Gives the first and second derivatives of the (cond) likelihood
//      w.r.t. the vector of betas.
//
//      Args:
//          rgivy (double *): Matrix of posteriors.
//			prior (double *): Matrix of priors.
//			x (double *): Matrix of indvars.
//          obs (int *): Number of observations.
//          classes (int *): Number of classes.
//          xcols (int *): Columns of x.
//      Return:
//          grad (double *): vector of first derivs.
//          hess (double *): Matrix of second derivs.
//
void d2lldbeta2_Ksites(int *numSite, double *rgivy, double *prior, double *x, int *obs,
                       int *classes, int *xcols, double *grad, double *hess) {
  int i,j,k,m,n,row,col,newrow,newcol,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int cxcols = *xcols;
  const int crank =  cxcols*(cclasses-1);
  const int csite = *numSite;
  
  for(s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      for (k=0;k<cxcols;k++) {
        //Gradient part
        for (j=0;j<(cclasses-1);j++) {
          grad[j*cxcols+k] += (double) x[k]*(rgivy[j]-prior[j]);
        }
        //Hessian part
        for (m=0;m<=k;m++) {
          for (j=0;j<(cclasses-1);j++) {
            col = j*cxcols + m;
            row = j*cxcols + k;
            //Diagonal block elements of the hessian
            hess[row*crank+col] += x[m]*x[k]*
              ( -prior[j]*(1.0-prior[j]) ) ;
            //In poLCA:
            //hess[row*crank+col] += x[m]*x[k]*
            //( -rgivy[j]*(1.0-rgivy[j]) +
            //prior[j]*(1.0-prior[j]) ) ;
            for (n=0;n<j;n++) {
              col = n*cxcols + m;
              //Subdiagonal elements of the hessian
              hess[row*crank+col] += x[m]*x[k]*
               ( prior[j]*prior[n] );
              //In poLCA:
              //hess[row*crank+col] += x[m]*x[k]*
                //( rgivy[j]*rgivy[n] - prior[j]*prior[n] );
            }
          }
        }
      }
      prior += cclasses;
      rgivy += cclasses;
      x += cxcols;
    }
  }
  
  // Copy the upper elements of the symetric off-diag blocks.
  for (i=0;i<(cclasses-1);i++) {
    for (j=i+1;j<(cclasses-1);j++) {
      for (k=0;k<cxcols;k++) {
        for (n=k+1;n<cxcols;n++) {
          row = j*cxcols+n;
          col = i*cxcols+k;
          newrow = j*cxcols+k;
          newcol = i*cxcols+n;
          hess[newrow*crank+newcol]=hess[row*crank+col];
        }
      }
    }
  }
  
  // Copy lower to upper off diagional elements
  for (col=0;col<(cclasses-1)*cxcols;col++) {
    for (row=0;row<col;row++) {
      hess[row*crank+col] = hess[col*crank+row];
    }
  }
}
//
void d2lldalpha2_Ksites(int *numSite, double *rgivy, double *prior, int *obs,
int *classes, double *grad, double *hess) {
  int i,j,m,row,col,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int csite = *numSite;
  const int crank =  csite*(cclasses-1);
  
  for (s=0;s<csite;s++){
  for (i=0;i<cobs;i++) {
      
      //Gradient part
      //grad: S * (C-1) length
      for (j=0;j<(cclasses-1);j++) {
        grad[j*csite+s] += (double) (rgivy[j]-prior[j]);
        
        // //Hessian part
        row = (cclasses-1)*s + j;
        //Subdiagonal elements of the hessian
        for (m=0;m<j;m++) {
          col = (cclasses-1)*s + m;
          hess[row*crank+col] += prior[j]*prior[m];
        }
        //Diagonal block elements of the hessian
        hess[row*crank+row] += -prior[j]*(1.0-prior[j]);
      }
    // Copy the upper elements of the symetric off-diag blocks.
    for (j=0;j<(cclasses-1);j++) {
      for (m=0;m<j;m++) {
        col = (cclasses-1)*s + j;
        row = (cclasses-1)*s + m;
        hess[row*crank+col]=hess[col*crank+row];
      }
    }
    prior += cclasses;
    rgivy += cclasses;
  }
}
}
//
void d2lldab2_Ksites(int *numSite, double *rgivy, double *prior, double *x, int *obs,
                     int *classes, int *xcols, double *hess) {
  int i,j,k,m,row,col,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int cxcols = *xcols;
  const int csite = *numSite;
  const int cranka = csite*(cclasses-1);
  
  for(s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      for(j=0;j<(cclasses-1);j++){
        for(k=0;k<cxcols;k++){
          row=j*cxcols+k;
          //alpha and beta in the same class
          hess[row*cranka+(s*(cclasses-1))+j] += -x[k]*(1-prior[j])*prior[j];
          // alapha and beta in different class
          for(m=0;m<j;m++){
            col=s*(cclasses-1)+m;
            hess[row*cranka+col] += x[k]*prior[j]*prior[m];
          }
          for(m=j+1;m<(cclasses-1);m++){
            col=s*(cclasses-1)+m;
            hess[row*cranka+col] += x[k]*prior[j]*prior[m];
          }
        }
      }
      prior += cclasses;
      rgivy += cclasses;
      x += cxcols;
    }
  }
}
//In Bandeen-Roche's paper, regression parameters beta are updated using L not Q
//Same gradient but different hessian
void d2lldbeta2_Ksites_BR(int *numSite, double *rgivy, double *prior, double *x, int *obs,
                       int *classes, int *xcols, double *grad, double *hess) {
  int i,j,k,m,n,row,col,newrow,newcol,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int cxcols = *xcols;
  const int crank =  cxcols*(cclasses-1);
  const int csite = *numSite;
  
  for(s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      for (k=0;k<cxcols;k++) {
        //Gradient part
        for (j=0;j<(cclasses-1);j++) {
          grad[j*cxcols+k] += (double) x[k]*(rgivy[j]-prior[j]);
        }
        //Hessian part
        for (m=0;m<=k;m++) {
          for (j=0;j<(cclasses-1);j++) {
            col = j*cxcols + m;
            row = j*cxcols + k;
            //Diagonal block elements of the hessian
            //In poLCA:
            hess[row*crank+col] += x[m]*x[k]*
            ( rgivy[j]*(1.0-rgivy[j]) -
            prior[j]*(1.0-prior[j]) ) ;
            for (n=0;n<j;n++) {
              col = n*cxcols + m;
              //Subdiagonal elements of the hessian
              //In poLCA:
              hess[row*crank+col] += x[m]*x[k]*
              ( -rgivy[j]*rgivy[n] + prior[j]*prior[n] );
            }
          }
        }
      }
      prior += cclasses;
      rgivy += cclasses;
      x += cxcols;
    }
  }
  
  // Copy the upper elements of the symetric off-diag blocks.
  for (i=0;i<(cclasses-1);i++) {
    for (j=i+1;j<(cclasses-1);j++) {
      for (k=0;k<cxcols;k++) {
        for (n=k+1;n<cxcols;n++) {
          row = j*cxcols+n;
          col = i*cxcols+k;
          newrow = j*cxcols+k;
          newcol = i*cxcols+n;
          hess[newrow*crank+newcol]=hess[row*crank+col];
        }
      }
    }
  }
  
  // Copy lower to upper off diagional elements
  for (col=0;col<(cclasses-1)*cxcols;col++) {
    for (row=0;row<col;row++) {
      hess[row*crank+col] = hess[col*crank+row];
    }
  }
}
//
void d2lldalpha2_Ksites_BR(int *numSite, double *rgivy, double *prior, int *obs,
                        int *classes, double *grad, double *hess) {
  int i,j,m,row,col,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int csite = *numSite;
  const int crank =  csite*(cclasses-1);
  
  for (s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      
      //Gradient part
      //grad: S * (C-1) length
      for (j=0;j<(cclasses-1);j++) {
        grad[j*csite+s] += (double) (rgivy[j]-prior[j]);
        
        // //Hessian part
        row = (cclasses-1)*s + j;
        //Subdiagonal elements of the hessian
        for (m=0;m<j;m++) {
          col = (cclasses-1)*s + m;
          hess[row*crank+col] += prior[j]*prior[m]-rgivy[j]*rgivy[m];
        }
        //Diagonal block elements of the hessian
        hess[row*crank+row] += -prior[j]*(1.0-prior[j])+rgivy[j]*(1.0-rgivy[j]);
      }
      // Copy the upper elements of the symetric off-diag blocks.
      for (j=0;j<(cclasses-1);j++) {
        for (m=0;m<j;m++) {
          col = (cclasses-1)*s + j;
          row = (cclasses-1)*s + m;
          hess[row*crank+col]=hess[col*crank+row];
        }
      }
      prior += cclasses;
      rgivy += cclasses;
    }
  }
}
//
void d2lldab2_Ksites_BR(int *numSite, double *rgivy, double *prior, double *x, int *obs,
                     int *classes, int *xcols, double *hess) {
  int i,j,k,m,row,col,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int cxcols = *xcols;
  const int csite = *numSite;
  const int cranka = csite*(cclasses-1);
  
  for(s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      for(j=0;j<(cclasses-1);j++){
        for(k=0;k<cxcols;k++){
          row=j*cxcols+k;
          //alpha and beta in the same class
          hess[row*cranka+(s*(cclasses-1))+j] += x[k]*((1.0-rgivy[j])*rgivy[j]-(1.0-prior[j])*prior[j]);
          // alapha and beta in different class
          for(m=0;m<j;m++){
            col=s*(cclasses-1)+m;
            hess[row*cranka+col] += x[k]*(prior[j]*prior[m]-rgivy[j]*rgivy[m]);  
          }
          for(m=j+1;m<(cclasses-1);m++){
            col=s*(cclasses-1)+m;
            hess[row*cranka+col] += x[k]*(prior[j]*prior[m]-rgivy[j]*rgivy[m]);  
          }
        }
      }
      prior += cclasses;
      rgivy += cclasses;
      x += cxcols;
    }
  }
}
//// Density Ratio tilting method
// tilt: T_k1
// postl : T_k1(X_1i) *w_kic(X_1i, Y_1i) l:leading
// probs: p^t : q first than class length=q*C
// probhat: class first than q
// For binary only!!!
void probhat_drt_Ksites(int *y, int *numSite, 
                        double *probs,double *postl, double *post,
                    int *items, int *obs, int *numChoices,
                    int *classes, double *probhat, double *denom) {
  //double *denom;
  int i,j,k,s;
  const int citems = *items;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int csite = *numSite;
  const double *firstprobs = probs;
  
  for (i=0;i<(citems*cclasses);i++) probhat[i]=0.0;
  
  denom = (double *) calloc((cclasses*citems),sizeof(double));
  for (i=0;i<(cclasses*citems);i++) denom[i]=0.0;
  
  for (s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      probs = (double *)firstprobs;
      for (j=0;j<cclasses;j++) {
        for (k=0;k<citems;k++) {
          if (y[k]>0) {
            probhat[j+cclasses*k] 
            += post[j]*(y[k]-1-probs[k]) + postl[j]*probs[k];
            denom[j*citems+k] += postl[j];
          }
        }
        probs += citems;
      }
      y+=citems;
      post+=cclasses;
      postl+=cclasses;
    }
  }
  
  for (j=0;j<cclasses;j++) {
    for (k=0;k<citems;k++) {
        probhat[j+cclasses*k] =
          (double) probhat[j+cclasses*k]/
            denom[j*citems+k];
    }
  }
  free(denom);
}
//tilt: T_k1(X_ki)
//priorl: lambda(X_1i)
//x: x1_all
void d2lldbeta2_drt_Ksites(int *numSite, double *tilt, double *rgivy, double *prior, 
                           double *priorl, double *x, double *xl,int *obs,
                       int *classes, int *xcols, double *grad, double *hess) {
  int i,j,k,m,n,row,col,newrow,newcol,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int cxcols = *xcols;
  const int crank =  cxcols*(cclasses-1);
  const int csite = *numSite;
  
  for(s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      //index = s*csite+i;
      for (k=0;k<cxcols;k++) {
        //Gradient part
        for (j=0;j<(cclasses-1);j++) {
          grad[j*cxcols+k] += (double) x[k]*(rgivy[j]-prior[j]);
        }
        //Hessian part
        for (m=0;m<=k;m++) {
          for (j=0;j<(cclasses-1);j++) {
            col = j*cxcols + m;
            row = j*cxcols + k;
            //Diagonal block elements of the hessian
            hess[row*crank+col] += tilt[0] * xl[m]*xl[k]*
              ( -priorl[j]*(1.0-priorl[j]) );
            for (n=0;n<j;n++) {
              col = n*cxcols + m;
              //Subdiagonal elements of the hessian
              hess[row*crank+col] += tilt[0] * xl[m]*xl[k]*
                ( priorl[j]*priorl[n] );
            }
          }
        }
      }
      prior += cclasses;
      priorl += cclasses;
      rgivy += cclasses;
      x += cxcols;
      xl += cxcols;
      tilt += 1;
    }
  }
  
  // Copy the upper elements of the symetric off-diag blocks.
  for (i=0;i<(cclasses-1);i++) {
    for (j=i+1;j<(cclasses-1);j++) {
      for (k=0;k<cxcols;k++) {
        for (n=k+1;n<cxcols;n++) {
          row = j*cxcols+n;
          col = i*cxcols+k;
          newrow = j*cxcols+k;
          newcol = i*cxcols+n;
          hess[newrow*crank+newcol]=hess[row*crank+col];
        }
      }
    }
  }
  
  // Copy lower to upper off diagional elements
  for (col=0;col<(cclasses-1)*cxcols;col++) {
    for (row=0;row<col;row++) {
      hess[row*crank+col] = hess[col*crank+row];
    }
  }
}
//
void d2lldalpha2_drt_Ksites(int *numSite, double *tilt, double *rgivy, double *prior, 
                            double *priorl, int *obs,
                        int *classes, double *grad, double *hess) {
  int i,j,m,row,col,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int csite = *numSite;
  const int crank =  csite*(cclasses-1);
  
  for (s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      
      //Gradient part
      //grad: S * (C-1) length
      for (j=0;j<(cclasses-1);j++) {
        grad[j*csite+s] += (double) (rgivy[j]-prior[j]);
        
        // //Hessian part
        row = (cclasses-1)*s + j;
        //Subdiagonal elements of the hessian
        for (m=0;m<j;m++) {
          col = (cclasses-1)*s + m;
          hess[row*crank+col] += tilt[0]*priorl[j]*priorl[m];
        }
        //Diagonal block elements of the hessian
        hess[row*crank+row] += -tilt[0]*priorl[j]*(1.0-priorl[j]);
      }
      // Copy the upper elements of the symetric off-diag blocks.
      for (j=0;j<(cclasses-1);j++) {
        for (m=0;m<j;m++) {
          col = (cclasses-1)*s + j;
          row = (cclasses-1)*s + m;
          hess[row*crank+col]=hess[col*crank+row];
        }
      }
      prior += cclasses;
      priorl += cclasses;
      rgivy += cclasses;
      tilt +=1;
    }
  }
}
void d2lldab2_drt_Ksites(int *numSite, double *tilt,
                         double *priorl, double *xl, int *obs,
                     int *classes, int *xcols, double *hess) {
  int i,j,k,m,row,col,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int cxcols = *xcols;
  const int csite = *numSite;
  const int cranka = csite*(cclasses-1);
  
  for(s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      for(j=0;j<(cclasses-1);j++){
        for(k=0;k<cxcols;k++){
          row=j*cxcols+k;
          //alpha and beta in the same class
          hess[row*cranka+(s*(cclasses-1))+j] += -tilt[1]*xl[k]*(1-priorl[j])*priorl[j];
          // alapha and beta in different class
          for(m=0;m<j;m++){
            col=s*(cclasses-1)+m;
            hess[row*cranka+col] += tilt[1]*xl[k]*priorl[j]*priorl[m];
          }
          for(m=j+1;m<(cclasses-1);m++){
            col=s*(cclasses-1)+m;
            hess[row*cranka+col] += tilt[1]*xl[k]*priorl[j]*priorl[m];
          }
        }
      }
      xl += cxcols;
      tilt +=1;
    }
  }
}
// Function: lambdahat
//     A function to return updates estimates the mixed proportion
//     for each class given the data and posterior distributions
//     Just a bunch of conditional means.
//
//  Args:
//     post (double *): Matrix of posterior estimates
//     obs (int *): Number of observations
//     classes (int):      Number of latent classes
void lambdahat_Ksites(int *numSite, double *post,
                      int *obs, int *classes, double *lambdahat) {
  int i,j,s;
  const int cobs = *obs;
  const int cclasses = *classes;
  const int csite = *numSite;
  
  for (s=0;s<csite;s++){
    for (i=0;i<cobs;i++) {
      for (j=0;j<cclasses;j++) {
        lambdahat[s*cclasses+j] += post[j]/cobs;
      }
      post+=cclasses;
    }
  }
}