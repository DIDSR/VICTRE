//#*********************************************************************************************************************************
//#*********************************************************************************************************************************
//#
//#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
//#								DIGITAL BREAST TOMOSYNTHESIS RECONSTRUCTION (FBP)
//#
//# BASED ON THE C RECONSTRUCTION CODE DEVELOPED BY LEESER, MUKHERJEE AND BROCK (BMC RESEARCH NOTES 7.1, p. 582, 2014).
//# WHICH IN TURN IMPLEMENTS THE FDK RECONSTRUCTION ALGORITHM DEVELOPED BY FESSLER (http://web.eecs.umich.edu/~fessler/).
//#
//# CODE MODIFIED BY: AUNNASHA SENGUPTA
//# 
//#
//# GITHUB LINK: 	https://github.com/DIDSR/VICTRE
//#
//#*********************************************************************************************************************************
//#*********************************************************************************************************************************

//##################################################################################################################################
//#
//#							DISCLAIMER
//#
//# This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. 
//# Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. 
//# Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, 
//# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, 
//# and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, 
//# its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. 
//# Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. 
//# Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, 
//# and any modified versions bear some notice that they have been modified. 
//#
//##################################################################################################################################



/*
 This program takes a number of 2D DBT projections and compute a 3D volume from those.
 It has four steps : Interpolation, Weighting, Filtering and Back projection.
 Usage: running the makefile should work
 */


#include "FBP_DBTrecon.h"
///////////////////////////////////////////////////////////////
//
// Cone beam CT back projection algorithm
//
// For details refer to Fessler's irt toolbox: cbct_back_mat.m
//
///////////////////////////////////////////////////////////////

void cbct_back(double * proj, int ns, int nt, int na, double * betas_rad, double * offset_xyz, int nx, int ny, int nz, double dx, double dy, double dz, double offset_s, double offset_t, double * img_t, int ia_skip, double dso, double dfs, double dsd, double ds, double dt, double orbit) {
   int idx, jdx, kdx = 0;
   double offset_source = 0;
   double * source_zs = (double*)malloc(sizeof(double)*na);
   memset(source_zs, 0, sizeof(double)*na);
   int scale_dang = 1;
   double orb; 
   clock_t tic, toc;
   double tspan = 0.0;


   double wx = (double(nx-1)/2.f) + offset_xyz[0];
   double wy = (double(ny-1)/2.f) + offset_xyz[1];
   double wz = (double(nz-1)/2.f) + offset_xyz[2];

   double * xc = (double*)malloc(sizeof(double)*nx*ny);
   double * yc = (double*)malloc(sizeof(double)*ny*nx-1);


   for (idx = 0; idx < ny; idx++) {
      int col = idx*nx;
      for (jdx = 0; jdx < nx; jdx++) {
      xc[col+jdx] = (double(jdx) - wx) * dx;
      }
   } 
   for (idx = 0; idx < ny; idx++) {
      int col = idx*nx;
      for (jdx = 0; jdx < nx; jdx++) {
    yc[col+jdx] = (double(idx) - wy) * dy;
      }
   }

   
   orb = orbit + fabs(betas_rad[1] - betas_rad[0]) ;
   double * zc = (double*)malloc(sizeof(double)*nz);
   for (idx = 0; idx < nz; idx++) {
      zc[idx] = (double(idx) - wz) * dz;
   }
 
   double ws = (double(ns+1)/2.f) + offset_s;
   double wt = (double(nt+1)/2.f) + offset_t;

   int sdim[2] = {0, 0};
   sdim[0] = ns+1; sdim[1] = nt;
   double * proj1 = (double*)malloc(sizeof(double)*sdim[0]*sdim[1]);
   memset(proj1, 0, sizeof(double)*sdim[0]*sdim[1]);

   for (idx = 0; idx < nz; idx++) {
   	int ia_min = 0;
        int ia_max = na;
        printf("\tImage slice %d ... \n", idx);
	tic = clock();

	// % loop over each projection angle
	double * img2 = (double*)malloc(sizeof(double)*nx*ny);
	memset(img2, 0, sizeof(double)*nx*ny);
        for (jdx = ia_min; jdx < ia_max; jdx += ia_skip) {
          double beta = betas_rad[jdx]*(M_PI/180.0);
          double * x_betas = (double*)malloc(sizeof(double)*nx*ny);
          for (kdx = 0; kdx < nx*ny; kdx++) {
	       x_betas[kdx] = (xc[kdx]*cosf(beta)) + (yc[kdx]*sinf(beta));
	       }
         double * y_betas = (double*)malloc(sizeof(double)*nx*ny);
         for (kdx = 0; kdx < nx*ny; kdx++) {
            y_betas[kdx] = dso - ((-1.f*xc[kdx]*sinf(beta)) + (yc[kdx]*cosf(beta)));
         }

          // % detector indices
          double * mag = (double*)malloc(sizeof(double)*nx*ny);
	  if ((dsd==INFINITY) || (dso==INFINITY)) {
          for (kdx = 0; kdx < nx*ny; kdx++) {
               mag[kdx] = 1.f;
               }
	     } 
          else {
	       for (kdx = 0; kdx < nx*ny; kdx++) {
		    mag[kdx] = dsd / y_betas[kdx];  
		    }
	       } 

          double * sprime = (double*)malloc(sizeof(double)*nx*ny);
          double * r_loop = (double*)malloc(sizeof(double)*nx*ny);
          if ((dfs==INFINITY) || (dsd==INFINITY) || (dso==INFINITY)) {
	  for (kdx = 0; kdx < nx*ny; kdx++) {
	  sprime[kdx] = mag[kdx] * x_betas[kdx];
	       }
	       } else if (dfs == 0.f) {
			 for (kdx = 0; kdx < nx*ny; kdx++) {
			      r_loop[kdx] = x_betas[kdx] - offset_source;
			      sprime[kdx] = dsd * atan2(r_loop[kdx], y_betas[kdx]);
		      }
		   }
		   double * tprime = (double*)malloc(sizeof(double)*nx*ny);
		   for (kdx = 0; kdx < nx*ny; kdx++) {
		      tprime[kdx] = mag[kdx] * (zc[idx] - source_zs[jdx]);
		   }
		   
         double * bs = (double*)malloc(sizeof(double)*nx*ny);
         for (kdx = 0; kdx < nx*ny; kdx++) {
             double tmp = sprime[kdx] / ds;
             bs[kdx] = tmp + ws;	   
             }
 

         double * bt = (double*)malloc(sizeof(double)*nx*ny);
         for (kdx = 0; kdx < nx*ny; kdx++) {
             double tmp = tprime[kdx] / dt;
             bt[kdx] = tmp + wt;
             }

         for (kdx = 0; kdx < nx*ny; kdx++) {
         if (bs[kdx] < 1 | bs[kdx] > ns) {bs[kdx]=ns+1;}
            }

         for (kdx = 0; kdx < nx*ny; kdx++) {
         if (bt[kdx] < 1){bt[kdx]=1;}
         else if (bt[kdx] > nt){bt[kdx]= nt;}
         }

        // % bi-linear interpolation:
	int * is = (int*)malloc(sizeof(int)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     is[kdx] = (int)floor(bs[kdx]);
	     }
	int * it = (int*)malloc(sizeof(int)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     it[kdx] = (int)floor(bt[kdx]);
	    }

	for (kdx = 0; kdx < nx*ny; kdx++) {
        if (is[kdx] == ns+1 ) {is[kdx]=ns;}
        if (it[kdx] == nt) {it[kdx]=nt-1;}
       	}

	// left weight
	double * wr = (double*)malloc(sizeof(double)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     wr[kdx] = bs[kdx] - (double)is[kdx];
	     }
	// right weight
	double * wl = (double*)malloc(sizeof(double)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     wl[kdx] = 1.f - wr[kdx];
	}
	// upper weight
	double * wu = (double*)malloc(sizeof(double)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     wu[kdx] = bt[kdx] - (double)it[kdx];
	     }
	// lower weight
	double * wd = (double*)malloc(sizeof(double)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     wd[kdx] = 1.f - wu[kdx];
	    }

	unsigned int * ibad = (unsigned int*)malloc(sizeof(unsigned int)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     ibad[kdx] = (is[kdx] < 0) | (is[kdx] > ns) | (it[kdx] < 0) | (it[kdx] > nt);
	     }
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     if (ibad[kdx]) {printf("BAD\n"); is[kdx] = ns+1; it[kdx] = nt+1; } 
	}


	int mdx = 0;
	int plane = ns*nt*jdx;
        for (kdx = 0; kdx < nt; kdx++) {
            int col = kdx*ns;
            int colx = (kdx)*sdim[0];
            for (mdx = 0; mdx < ns; mdx++) {
                proj1[colx+mdx] = proj[plane+col+mdx];
            }
         }
      
         
	// vertical interpolation
	double * p1 = (double*)malloc(sizeof(double)*nx*ny);
	double * p2 = (double*)malloc(sizeof(double)*nx*ny);
	double * p0 = (double*)malloc(sizeof(double)*nx*ny);
	for (kdx = 0; kdx < nx*ny; kdx++) {
	     int is1 = is[kdx]-1; int is2 = is1 + 1;
             int it1 = it[kdx]-1; int it2 = it1 + 1;
	     p1[kdx] = (wl[kdx] * proj1[((it1)*sdim[0])+is1]) + (wr[kdx] * proj1[((it1)*sdim[0])+is2]);
	     p2[kdx] = (wl[kdx] * proj1[((it2)*sdim[0])+is1]) + (wr[kdx] * proj1[((it2)*sdim[0])+is2]);
	     p0[kdx] = (wd[kdx] * p1[kdx]) + (wu[kdx] * p2[kdx]);                       // corrected equation
		}
	if ((dfs == INFINITY) || (dsd == INFINITY) || (dso == INFINITY)) {
	    for (kdx = 0; kdx < nx*ny; kdx++) {
	         p0[kdx] = p0[kdx] * mag[kdx] * mag[kdx];
                 }
	     } else if (dfs == 0.f) {
            for (kdx = 0; kdx < nx*ny; kdx++) {
		 p0[kdx] = (p0[kdx] * (dsd*dsd)) / ((r_loop[kdx]*r_loop[kdx]) + (y_betas[kdx]*y_betas[kdx]));
	        }
         }

        for(kdx = 0; kdx < nx*ny; kdx++) {
            img2[kdx] = img2[kdx] + p0[kdx];
	    }
	   
	 free(x_betas); free(y_betas); free(mag);
         free(sprime); free(tprime); free(r_loop);
         free(bs); free(bt);
         free(is); free(it);
         free(wr); free(wl); free(wu); free(wd);
         free(ibad);   
         free(p1); free(p2); free(p0);
      } // end ia (projections angles)
      
      int plane = nx*ny*idx;
      for (kdx = 0; kdx < nx*ny; kdx++) {        
        img_t[plane+kdx] = img2[kdx];
      }
      free(img2);
      
      if (scale_dang) {
      for (kdx = 0; kdx < nx*ny; kdx++) { img_t[plane+kdx] =  img_t[plane+kdx] * ( (0.5*(fabs(orb)*M_PI/180.0)) / (double(na)/double(ia_skip)) ); }
      } 
      toc = clock();
      tspan = (double)(toc-tic) / (double)CLOCKS_PER_SEC;
   } // end iz (slices)

   free(source_zs);
   free(xc); free(yc); free(zc);
   free(proj1);
   return;
}
/////////////////////////////////////////For details refer to Fessler's irt toolbox: feldkamp_weight1.m //////////////////////
//
// Feldkamp Weight algorithm
//
void fdk_weight(  double * proj, int ns, int nt, int na, double ds, double dt, double offset_s,
                  double offset_t, double dsd, double dso, double dfs, int w1cyl) {
   int idx, jdx, kdx, offset, index = 0;
   double tmp1, tmp2, tmp3 = 0.f;
   double * ss = (double*)malloc(sizeof(double)*ns*nt);

   for (idx = 0; idx < ns; idx++) {
      tmp1 = ((idx - ((double(ns)-1.f) / 2.f)) - offset_s) * ds;
      for (jdx = 0; jdx < nt; jdx++) { ss[(jdx*ns)+idx] = tmp1; }
   }
   double * tt = (double*)malloc(sizeof(double)*ns*nt);
   for (idx = 0; idx < nt; idx++) {
      tmp2 = ((idx - ((double(nt)-1.f) / 2.f)) - offset_t) * dt;
      offset = idx*ns;
      for (jdx = 0; jdx < ns; jdx++) { tt[offset+jdx] = tmp2; }
   }
   
   double * ww = (double*)malloc(sizeof(double)*ns*nt);
   if (dfs == INFINITY) {
      if (w1cyl){   
      for (idx = 0; idx < nt; idx++) {
         offset = idx*ns;
         for (jdx = 0; jdx < ns; jdx++) {
            index = offset+jdx;
            tmp1 = ss[index]*ss[index];
            tmp2 = tt[index]*tt[index];
	    tmp3 = (tt[index]/dsd)*(tt[index]/dsd);
            ww[offset+jdx] = (dso / sqrt(dsd*dsd + tmp1 + tmp2));
         }
       }
      }
   /////////more added conditions for when the weighting is not exact for non-cylindrical objects //////////////////
      else if (!w1cyl) {
         for (idx = 0; idx < nt; idx++) {
            offset = idx*ns;
            for (jdx = 0; jdx < ns; jdx++) {
               index = offset + jdx;
               tmp1 = (tt[index] / dsd);
               ww[index] *= sqrt(1.f + (tmp1*tmp1));
            }
         }
      }
	} else if (dfs == 0.f) {
           if (!w1cyl){
	   for (idx = 0; idx < nt; idx++) {
	      offset = idx*ns;
	      for (jdx = 0; jdx < ns; jdx++) {
	         index = offset + jdx;
	         tmp1 = (tt[index] / dsd);
	         ww[index] = (dso/dsd) * cosl(ss[index] / (dsd * sqrt(1.f + (tmp1*tmp1))));
	      }
	   }
	}
	   else if (w1cyl) {
	      for (idx = 0; idx < nt; idx++) {
	         offset = idx*ns;
	         for (jdx = 0; jdx < ns; jdx++) {
	            index = offset + jdx;
	            tmp1 = (tt[index] / dsd);                                       
                    ww[index] = (dso/dsd) * cosl(ss[index]/dsd)/sqrt(1.f + (tmp1*tmp1));           // corrected equation
                  
	         }
	      }
	   }
	} else {
	   printf("ERROR: Other configurations not implemented\n");
   }
   for (idx = 0; idx < na; idx++) {
      offset = idx*ns*nt;
      for (jdx = 0; jdx < nt; jdx++) {
         index = offset + (jdx*ns);
         for (kdx = 0; kdx < ns; kdx++) {
            proj[index+kdx] *= ww[(jdx*ns)+kdx];
         }
      }
   }

   free(ss);
   free(tt);
   free(ww);
   return;
}

//
// Ramp Arc function
//
void ramp_arc(double * hval, int n, double ds, double dsd) {
   int idx = 0;
   double tmp = 0;
   int * nn = (int*)malloc(sizeof(int)*n);
   for (idx = 0; idx < n; idx++) {
      nn[idx] = idx - (n/2);
      if (nn[idx] == 0) {
         hval[idx] = 1.f / (4.f*ds*ds);
      } else {
         hval[idx] = 0.f;
      }
   }
   for (idx = 0; idx < n; idx++) {
      if ((nn[idx]%2) != 0) {
         tmp = sin(((double)(nn[idx]) * ds) / dsd) * dsd * M_PI;
         hval[idx] = -1.f / (tmp*tmp);
      }
   }
   free(nn);
   return;
}

//
// Ramp filter function
//
void ramp_flat(double * hval, int n, double ds) {
   int idx = 0;
   int * nn = (int*)malloc(sizeof(int)*n);
   for (idx = 0; idx < n; idx++) {
      nn[idx] = idx - (n/2);
      hval[idx] = 0.0;
   }
   hval[(n/2)] = 1.0 / 4.0;
   for (idx = 0; idx < n; idx++) {
      if ((abs(nn[idx])%2) == 1) {
         hval[idx] = -1.0/ ((double)(M_PI * nn[idx]) * (M_PI * nn[idx]));
         
      }
        hval[idx] = hval[idx]/(ds*ds);
   }
   free(nn);
   return;
}

void fbp_ramp(double * hval, char * type, int n, double ds, double dsd) {

	if (!strcmp(type,"arc")) {
	   ramp_arc(hval, n, ds, dsd);
   } else if (!strcmp(type,"flt")) {
	   ramp_flat(hval, n, ds);
	} else {
	   printf("ERROR: bad fan type\n");

   }
   return;
}

////////////////////////////////////For details refer to Fessler's irt toolbox: fdk_fan_filter.m //////////////////////
//
// Feldkamp-Davis-Kress fan filter function
// 
void fdk_fan_filter(double * hval, char * type, int n, double ds, double dsd, char * wintype ) {
   int idx = 0;
   int offset = n/2;
   float * h_shift = (float*)malloc(sizeof(float)*n);
   int w;
   double ii;
   fftwf_complex * fftc_out = (fftwf_complex*)fftwf_malloc(sizeof(fftw_complex)*((n/2)+1));
   fftwf_plan plan = fftwf_plan_dft_r2c_1d(n, h_shift, fftc_out, 0);

   fbp_ramp(hval, type, n, ds, dsd);

   for (idx = 0; idx < n; idx++) {
      if (idx < (n/2)) {
         h_shift[idx] = hval[idx+offset];
      } else {
         h_shift[idx] = hval[idx-offset];
      }
   }

   fftwf_execute(plan);

   for (idx = 0; idx < (n/2)+1; idx++) {
      hval[idx] = fftc_out[idx][0];
   }
   for (idx = 0; idx < (n/2); idx++) {
      hval[idx+(n/2)+1] = fftc_out[(n/2)-1-idx][0];
   }

  
   double * window = (double*)malloc(sizeof(double)*n);
	if (!strcmp(wintype, "ramp")) {
	   for (idx = 0; idx < n; idx++) { window[idx] = 1.f; }
   } else if (!strcmp(wintype,"hann")) {
      for (idx = 0; idx < n; idx++) { window[idx] = 0.5 * (1.0 - cosl(((double)2.0*M_PI*idx)/(n))); }
	}
     else if (!strcmp(wintype,"han7")) { w = round(n*0.75);                  //Filter added to reduce noise in reconstructed slices
      for (idx = 0; idx<n; idx++){ ii = idx-offset;
        if (fabs(ii)<(w/2)){
        window[idx] = 0.5*(1.0 + cosl(((double)2.0*M_PI*ii)/w));}
        else {window[idx] = 0;} 
      }
        }
     else {
      printf("ERROR: Unknown window type!\n");
   }


   for (idx = 0; idx < n; idx++) {
      if (idx < (n/2)) {
         h_shift[idx] = window[idx+offset];
      } else {
         h_shift[idx] = window[idx-offset];
      }
      hval[idx] *= h_shift[idx];
      
   }


   fftwf_destroy_plan(plan);
   fftwf_free(fftc_out);
   fftwf_free(h_shift);
   return;
}

/////////////////////////////////////////For details refer to Fessler's irt toolbox: fdk_filter.m //////////////////////
//
// Feldkamp-Davis-Kress filter function
//
void fdk_filter(  double * proj, char * window, double dsd, double dfs, double ds,
                  int ns, int nt, int na) {
   int idx, jdx, kdx, offset, index = 0;
   FILE * fid = NULL;
   int npad = 1<<(int)ceil(log2((2*ns)-1));
   float * proj_pad = (float*)malloc(sizeof(float)*npad);
   for (idx = 0; idx < npad; idx++) { proj_pad[idx] = 0.f; }
   fftwf_complex * fftout1 = (fftwf_complex*)fftwf_malloc(sizeof(fftw_complex)*(npad));
   fftwf_plan planF = fftwf_plan_dft_r2c_1d(npad, proj_pad,fftout1, 0);
   fftwf_plan planR = fftwf_plan_dft_c2r_1d(npad, fftout1, proj_pad, 0);

   double * hval  = (double*)malloc(sizeof(double)*npad);
   for (idx = 0; idx < npad; idx++) { hval[idx] = 0.f; }
   if (dsd == INFINITY) {
      fdk_fan_filter(hval, "flt", npad, ds, 0.f, window);
	} else if (dfs == INFINITY) {
      fdk_fan_filter(hval, "flt", npad, ds, 0.f, window);
	} else if (dfs == 0.f) {
      fdk_fan_filter(hval, "arc", npad, ds, dsd, window);
   }
   for (idx = 0; idx < npad; idx++) { hval[idx] *= ds; }

   for (idx = 0; idx < na; idx++) {
      offset = idx*ns*nt;
      for (jdx = 0; jdx < nt; jdx++) {
         index = offset+(jdx*ns);
         for (kdx = 0; kdx < npad; kdx++) { 
         if (kdx < ns)
         {proj_pad[kdx] = (float)proj[index+kdx];}
         else 
         { proj_pad[kdx] = 0.f; }
         }
          
         for (kdx = 0; kdx < npad; kdx++) { fftout1[kdx][0] = 0.f;fftout1[kdx][1] = 0.f; }
 
         fftwf_execute(planF);

         for (kdx = 0; kdx < npad; kdx++) {
            fftout1[kdx][0] *= hval[kdx];
            fftout1[kdx][1] *= hval[kdx];
         }
         fftwf_execute(planR);
         for (kdx = 0; kdx < npad; kdx++) { proj_pad[kdx] /=  (double)(npad); }
         for (kdx = 0; kdx < ns; kdx++) { proj[index+kdx] = proj_pad[kdx]; }
 
     }
   }

   fftwf_destroy_plan(planF);
   fftwf_destroy_plan(planR);
   fftwf_free(fftout1);
   free(proj_pad);
   free(hval);
   return;
}
void freeArray(double **a, int m) {
    int i;
    for (i = 0; i < m; ++i) {
        free(a[i]);
    }
    free(a);
}

void freeArrayi(int **a, int m) {
    int i;
    for (i = 0; i < m; ++i) {
        free(a[i]);
    }
    free(a);
}

/////////////////////For Details Refer to MATLAB code by Dr. Rongping Zeng : fbp_dbt.m////////////////////////////////
void ct_projections(  double * projv, double * proji, double * betas_rad, int i, int ns, int nt, int ns_old, int nt_old, double offset_s, double offset_t, double offset_sold, double offset_told, double dt, double ds, double dod, double dso, double dsd, double orbit){
	
	double halforbit, cbg_det_wid, cbg_det_depth, beta, F1, F2, F3, F4, ang, acos, asin, ws, wt, ws_old, wt_old, xs, ys, zs = 0;
        int idx, jdx, nrows, ncols, flagds, flagdt, dis, djs, dit, djt;
  	FILE * fid = NULL;
        beta = betas_rad[i];
        ang = (beta*M_PI)/180;
        acos = cos(ang);
        asin = sin(ang);
        ws = ((double)ns-1)/2 + offset_s;
        wt = ((double)nt-1)/2 + offset_t;

          

        ws_old = ((double)ns_old-1)/2 + offset_sold;
        wt_old = ((double)nt_old-1)/2 + offset_told;
        xs = -dso*asin;
        ys = 0;
        zs = dso*acos;
        ncols = ns_old;
	nrows = nt_old;
        
        double * sn = (double*)malloc(sizeof(double)*ns);
        memset(sn, 0, ns);
        double * tn = (double*)malloc(sizeof(double)*nt);
        memset(tn, 0, nt);

        
        double * s_old = (double*)malloc(sizeof(double)*ns_old);
        memset(s_old, 0, ns_old);
        double * t_old = (double*)malloc(sizeof(double)*nt_old);
        memset(t_old, 0, nt_old);
        double **xo = (double **)malloc(sizeof(double *)*nt_old);
        double **yo = (double **)malloc(sizeof(double *)*nt_old);
  
        
        for (idx=0; idx<nt_old; idx++)
          {
                xo[idx] = (double *)malloc(ns_old * sizeof(double));
                yo[idx] = (double *)malloc(ns_old * sizeof(double));
                memset(xo[idx], 0, ns_old);
                memset(yo[idx], 0, ns_old);
          }

        
        for (idx = 0; idx<ns; idx++)
           {
           	sn[idx] = ds*((double)(idx)-ws);
                
           }
        for (idx = 0; idx<nt; idx++)
           {
           	tn[idx] = dt*((double)(idx)-wt);
           }
        for (idx = 0; idx<ns_old; idx++)
           {
           	s_old[idx] = ds*((double)(idx)- ws_old);
           }
        for (idx = 0; idx<nt_old; idx++)
           {
           	t_old[idx] = dt*((double)(idx)- wt_old);
           }


        for (idx = 0; idx<nt_old; idx++){
    	   for (jdx = 0; jdx<ns_old; jdx++){
        	xo [idx][jdx] = s_old[jdx];
        	yo [idx][jdx] = t_old[idx];
	    }
          }
        
        free(s_old); free(t_old);
        
        double **k = (double **)malloc(sizeof(double *)*nt);
        double **xi = (double **)malloc(sizeof(double *)*nt);
        double **yi = (double **)malloc(sizeof(double *)*nt);
        
        for (idx=0; idx<nt; idx++)
          {
                k[idx]  = (double *)malloc(ns * sizeof(double));
                memset(k[idx], 0, ns);
                xi[idx] = (double *)malloc(ns * sizeof(double));
                memset(xi[idx], 0, ns);
                yi[idx] = (double *)malloc(ns * sizeof(double));
                memset(yi[idx], 0, ns);
          }
       
      for (idx = 0; idx<nt; idx++){
         for (jdx = 0; jdx<ns; jdx++){
        	k[idx][jdx]  = (-dod-(sn[jdx]*asin - dod*acos))/(zs-(sn[jdx]*asin - dod*acos));
        	xi[idx][jdx] = (sn[jdx]*acos + dod*asin) + k[idx][jdx]*(xs-(sn[jdx]*acos + dod*asin));
        	yi[idx][jdx] = tn[idx] + k[idx][jdx]*(ys-tn[idx]);
       }
     }
      freeArray(k, nt);
 

    	double **s = (double **)malloc(sizeof(double *)*nt);
        double **t = (double **)malloc(sizeof(double *)*nt);
       
        for (idx=0; idx<nt; idx++)
          {
                s[idx] = (double *)malloc(ns * sizeof(double));
                t[idx] = (double *)malloc(ns * sizeof(double));
                memset(s[idx], 0, ns);
                memset(t[idx], 0, ns);
          }
       for (idx = 0; idx<nt; idx++){
         for (jdx = 0; jdx<ns; jdx++){
             s[idx][jdx] = 1+ ((xi[idx][jdx]-xo[0][0])/(xo[nt_old-1][ns_old-1]-xo[0][0])*(ncols-1));
             t[idx][jdx] = 1+ ((yi[idx][jdx]-yo[0][0])/(yo[nt_old-1][ns_old-1]-yo[0][0])*(nrows-1));
       }
     }

        freeArray(xi, nt); freeArray(yi, nt); freeArray(xo, nt_old); freeArray(yo, nt_old);
        

        int **ndx = (int **)malloc(sizeof(int *)*nt);
        double **flags = (double **)malloc(sizeof(double *)*nt);
        double **flagt = (double **)malloc(sizeof(double *)*nt);
        
        for (idx=0; idx<nt; idx++)
          {
                ndx[idx] = (int *)malloc(ns * sizeof(int));
                flagt[idx] = (double *)malloc(ns * sizeof(double));
		flags[idx] = (double *)malloc(ns * sizeof(double));
          }
      

	for( idx =0;idx<nt; idx++){
    	    for (jdx = 0; jdx<ns; jdx++) {
        flagds = 0;
        flagdt = 0;
        flagt[idx][jdx]  = 0;
        flags[idx][jdx]  = 0;
        dis = 0;
        djs = 0;
	dit = 0;
        dis = 0; 
        if ((s[idx][jdx] < 1)||(s[idx][jdx] > (ncols))){
           s[idx][jdx] = 1;
           flags[idx][jdx] = 1;
        }
        if ((t[idx][jdx] < 1) || (t[idx][jdx] > (nrows)) ){
           t[idx][jdx] = 1;
           flagt[idx][jdx]= 1;
        }
        
        ndx[idx][jdx] = floor(t[idx][jdx])+floor(s[idx][jdx]-1)*nrows;
        
       // Compute intepolation parameters, check for boundary value. 
            if (s[idx][jdx] == (ncols)){
                dis = idx;
                djs = jdx;
                flagds = 1;
              }
        
        s[idx][jdx] = (s[idx][jdx] - floor(s[idx][jdx]));
        
        if (flagds == 1){
        s[dis][djs] = s[dis][djs] +1 ;
        ndx[dis][djs] = ndx[dis][djs] - nrows ;
        }
        
       if (t[idx][jdx] == (nrows)){
                dit = idx;
                djt = jdx;
                flagdt = 1;
       }
        
        t[idx][jdx] = (t[idx][jdx] - floor(t[idx][jdx]));

        if (flagdt == 1){
        t[dit][djt] = t[dit][djt] +1 ;
        ndx[dit][djt] = ndx[dit][djt] - 1 ;
        }
         }
        }
 	
        double **onemt = (double **)malloc(sizeof(double *)*nt);
        double **F = (double **)malloc(sizeof(double *)*nt);
       
        for (idx=0; idx<nt; idx++)
          {
                onemt[idx] = (double *)malloc(ns * sizeof(double));
                F[idx] = (double *)malloc(ns * sizeof(double));
                memset(onemt[idx], 0, ns);
                memset(F[idx], 0, ns);
          }

	for (idx = 0; idx<nt; idx++){
   	 for (jdx = 0;jdx<ns; jdx++){
        //Now interpolate.
        onemt[idx][jdx] = 1-t[idx][jdx];
        F1 = projv[ndx[idx][jdx]-1]*(onemt[idx][jdx]);
        F2 = projv[ndx[idx][jdx]]*t[idx][jdx] ;
        F3 = projv[ndx[idx][jdx]+nrows-1]*(onemt[idx][jdx]);
        F4 = projv[ndx[idx][jdx]+(nrows)]*t[idx][jdx];
        F[idx][jdx] =  ( F1 + F2)*(1-s[idx][jdx]) + ( F3 + F4 )*s[idx][jdx];
    
         if ((flags [idx][jdx] == 1) ||(flagt[idx][jdx] == 1) ){
            F[idx][jdx] = 0;
         }  
        }
       }
     
      
     for (idx = 0; idx<nt; idx++){
      for (jdx = 0; jdx<ns; jdx++){
        proji[idx*ns + jdx] = F[idx][jdx];
      }
     }

   freeArray(F, nt); 
   freeArray(onemt, nt); freeArray(s, nt); freeArray(t, nt); freeArray(flags, nt); freeArray(flagt, nt); freeArrayi(ndx, nt); 
  
return;
}


// 
// Main CBCT FUNCTION 
// 
int main(int argc, char * argv[]) 
{
  // CBCT Variables
  int ind, plane, index, indexi, ns_old, nt_old, na, nx, ny, nz, down, phntm  = 0;
  int i, j, x, z, a, b, c, di, dj, jdx, kdx, nt, ns, flat_corr;
  double zfov, orbit, ds, dt, dx, dy, dz, dso, dfs, dod, dsd = 0.f;
  int f;
  double d_objbottom_det, s_old_end, t_old_end, halforbit, cbg_det_wid, cbg_det_depth, ws_old, wt_old= 0;
  double phantom_voxsz, phantom_szy, phantom_szx, phantom_szz;
  int phantom_dimx, phantom_dimy, phantom_dimz;
  double offset_s, offset_t = 0.f;
  double offset_sold, offset_told = 0.f;
  int ia_skip, w1cyl, idx = 0;
  double offset_xyz[3] = {0.f, 0.f, 0.f};
  char window[5] = {'h', 'a', 'n', '7', '\0'};
  double *proji, *temp, *temp2, *proj, *projdbt, *betas_rad,*b_rad, *img, *img_temp,  *img_tr, *img_cmp = NULL;
  float *projdbt_un,*flatfield ;
  double *angles;
  FILE * fid = NULL;
  clock_t tic, tac, toc = 0;
  double tinc, tspan = 0.0;
  double orb_start, del_orb; 
  char inFile[MAXSTRLEN], ffFile[MAXSTRLEN];
  char outFile[MAXSTRLEN];
  char ph_seed[20];
  uint phant = 0;

  // Start timer
  tic = clock();



  /////////////////////////////////////////////// Reading Command Line Arguments //////////////////////////////////////////

  if(argc != 19)
    { printf("\n\n\t Missing arguments!! Program exiting!! \n\n"); exit(0);}

  na = atoi(argv[1]);
  ns_old = atoi(argv[2]);
  nt_old = atoi(argv[3]);
  ds = atof(argv[4]);
  dsd = atof(argv[5]);
  dso = atof(argv[6]);
  offset_sold = atof(argv[7]);
  orbit = atof(argv[8]);
  phantom_dimx = atoi(argv[9]);
  phantom_dimy = atoi(argv[10]);
  phantom_dimz = atoi(argv[11]);
  phantom_voxsz = atof(argv[12]);
  dx = atof(argv[13]);
  dz = atof(argv[14]);
  offset_xyz[0] = atof(argv[15]);
  orb_start = atof(argv[16]);
  del_orb = atof(argv[17]);
  strcpy(ph_seed, argv[18]); // random seed from VICTRE phantom file name                                                                                                                                       

  printf("\n Command line arguments:");
  printf ("\n\n na = %i \n ns = %i \n nt = %i \n ds = %lf \n dsd = %lf \n dso = %lf \n offset_s = %lf \n orbit = %lf \n phantom_dimX = %i \n phantom_dimY = %i \n phantom_dimZ = %i \n phantom_voxsz = %lf \n dx = %lf \n dz = %lf \n offset_x = %lf \n orbit start angle = %lf \n delta angle = %lf \n seed = %s \n", na, ns_old, nt_old, ds, dsd, dso, offset_sold, orbit, phantom_dimx, phantom_dimy, phantom_dimz, phantom_voxsz, dx, dz, offset_xyz[0], orb_start, del_orb, ph_seed);


  ia_skip = 1;                            
  w1cyl   = 1;                          
  flat_corr = 1;                      //flag to correct with flat field
  down = 1;                          //downsampling factor: for testing purposes down=4


/////////////////////////SCANNER GEOMETRY ////////////////////// 
  ns_old = ns_old/down;    //number of detector elements in the 's' direction                    
  nt_old = nt_old/down;    //number of detector elements in the 't' direction
  ds = ds*down;  //    in cm: detector element pixel size; 
               //'s' is the direction of the x-ray tube moving direction 
              // and 't' is the perpendicular direction to 's' direction. 
  dt = -ds;
  dod = dsd-dso;    //in cm: dist. from the rotation ctr to the detector
  dfs = INFINITY;   // Infinite for flat detector
  offset_told = -((double)nt_old)/2; //detector center offset along the 't' direction in pixels relative to the x-ray source center
  d_objbottom_det = 2-(0.1/2); //in cm, distance from the bottom of the object to the detector

///////////////RECON VOLUME GEOMETRY ////////////////////////////
// voxel size (dx, dy, dz) in cm, 
// dimensions (nx, ny, nz) and 
// center offsets (offset_x, offset_y, offset_z) in pixels relative to the rotation center.

  phantom_szx =  phantom_dimx*phantom_voxsz;
  phantom_szy =  phantom_dimy*phantom_voxsz; 
  phantom_szz =  phantom_dimz*phantom_voxsz;

  dx = dx*down;
  dy = dx;
  dz = dz*down;

  nx = ceil(phantom_szx/dx);
  ny = ceil(phantom_szy/dy);
  nz = ceil(phantom_szz/dz);
  printf("nx = %i\t ny = %i\t nz = %i\n", nx,ny,nz);
  offset_xyz[2] = -((double)ny)/2;    // 0 for full cone, -(ny)/2 for half cone
  zfov = nz*dz;

  offset_xyz[1] = (dod-((zfov/2)+d_objbottom_det))/dz;  //in pixels: offset of the volume ctr to the rotation ctr in the z direction; 
                                                       //downward direction is the positive z direction
                                                       //The provided value here is when the object was placed right above the detector.

///////////////FLAT FIELD FILE ////////////////////////////

  strcpy(ffFile,"PATH FOR THE CONCATENATED FLATFIELD PROJECTIONS RAW FILE/flatfield_25proj.raw");

///////////////PHANTOM DATA///////////////////////////////


  for (phntm = 0; phntm < 1; phntm++)      // run for phantom containing ph_seed
    {   
      //      fscanf (fp, "%s\t%s\n", inFile, outFile);
      sprintf(inFile, "PATH TO THE CONCATENATED DBT PROJECTIONS FILE/DBT_%s_25proj.raw", ph_seed);      //set path for the MC-GPU input projections
      sprintf(outFile, "PATH TO SAVE THE OUTPUT RECONSTRUCTED FILE/DBT_%s_recon.raw", ph_seed); //set path for the reconstructed output volumes
      printf("PHANTOM NAME %s\n", inFile);
      printf("O/P FILENAME NAME %s\n", outFile);

      // Read input variables from file based on run level
      
      printf("Reading input data ... \n");
      if (runlvl == PREWEIGHT) 
	{
	  fid = fopen(inFile, "rb");
	  if (fid == NULL) { printf("\nERROR: Could not open pre-weight file %s\n", inFile); }
	} 
      else if (runlvl == PREFILTER) 
	{
	  fid = fopen(inFile, "rb");
	  if (fid == NULL) { printf("\nERROR: Could not open pre-filter file %s\n", inFile); }
	} 
      else if (runlvl == BACKPROJ) 
	{
	  fid = fopen(inFile, "rb");
	  if (fid == NULL) { printf("\nERROR: Could not open back projection file %s\n", inFile); }
	}




/////////////////////////////////////////////////////////////////
  
      if (runlvl <= PREWEIGHT) 
	{
	  projdbt = (double*)malloc(sizeof(double)*ns_old*nt_old*na);
	  projdbt_un = (float*)malloc(sizeof(float)*ns_old*nt_old*na*down*down);
	  temp2 = (double*)malloc(sizeof(double)*ns_old*nt_old*na*down*down);
	  fread(projdbt_un, sizeof(double), ns_old*nt_old*na*down*down, fid);
	  betas_rad = (double*)malloc(sizeof(double)*na);
	  fclose(fid);  
	  betas_rad[0] = orb_start;
	  for (x=1; x<na; x++)
	    {
	      betas_rad[x]= betas_rad[x-1] + del_orb;
	    }
	  printf("DONE READING INPUT FOR THIS PHANTOM!\n"); 

	  // Flat field corections 
	  if ((flat_corr == 1))
	    {
	      fid = fopen(ffFile, "rb");
	      if (fid == NULL) { printf("\nERROR: Could not open flat field file %s\n", ffFile); }
	      flatfield = (float*)malloc(sizeof(float)*ns_old*nt_old*na*down*down);
	      fread(flatfield, sizeof(float), ns_old*nt_old*na*down*down, fid);
	      fclose(fid);
	      printf("DONE READING INPUTS!\n");
    
	      for (kdx = 0; kdx < nt_old*ns_old*na*down*down; kdx ++ )
		{
		  temp2[kdx] =  log (flatfield[kdx]) - log (projdbt_un[kdx]);
		  if (temp2[kdx]<0) {temp2[kdx] = 0;}
		  if ((isnan(temp2[kdx])) | (isinf(temp2[kdx]))) {temp2[kdx] = 0;}
		}
	      printf("DONE CORRECTING WITH FLAT FIELD!\n");
	      for (kdx = 0; kdx< na; kdx++)
		{
		  ind   = kdx*nt_old*ns_old*down*down; 
		  plane = ((na-kdx)-1)*nt_old*ns_old;  
		  for (idx=0; idx<nt_old*down; idx=idx+down)
		    {
		      for (jdx = 0; jdx<ns_old*down; jdx=jdx+down)
			{
			  projdbt[plane+(ns_old)*((idx)/down)+ ((jdx)/down)] = temp2[ind+ns_old*down*idx+ jdx]; 
			}
		    }
		}
	      
	    }

	} // runlvl <= PREWEIGHT if ends

      // Initialize output image and intermediate outputs
      img = (double*)malloc(sizeof(double)*nx*ny*nz);
      img_temp = (double*)malloc(sizeof(double)*nx*ny*nz);
      img_tr = (double*)malloc(sizeof(double)*nx*ny*nz);
      for (idx = 0; idx < nx*ny*nz; idx++) { img[idx] = 0.f; }
      for (idx = 0; idx < nx*ny*nz; idx++) { img_temp[idx] = 0.f; }
      for (idx = 0; idx < nx*ny*nz; idx++) { img_tr[idx] = 0.f; }

      //  Manipulating DBT into CT goemetry 
      ws_old = ((double)ns_old-1)/2 + offset_sold;
      wt_old = ((double)nt_old-1)/2 + offset_told;
      s_old_end = ds*((double)(ns_old-1)- ws_old);
      t_old_end = dt*((double)(nt_old-1)- wt_old);
      
      halforbit = (orbit/2)*(M_PI/180);
      cbg_det_wid = dod*tan(halforbit) + s_old_end/cos(halforbit);
      cbg_det_depth = t_old_end*(dsd+s_old_end*tan(halforbit))/dsd;
      ns = ceil(cbg_det_wid/ds)*2;
      nt = ceil(cbg_det_depth/(dt))*2;

      offset_t = 0.0;

      proji = (double*)malloc(sizeof(double)*ns*nt);
      proj = (double*)malloc(sizeof(double)*ns*nt*na);
    
      if (runlvl <= PREWEIGHT) 
	{
	  temp = (double*)malloc(sizeof(double)*ns_old*nt_old);
	  for (idx = 0; idx < na*ns*nt; idx++) { proj[idx]= 0.f; }
	  for (kdx = 0; kdx< na; kdx++)
	    {
	      index = kdx*nt_old*ns_old;   
	      for (idx=0; idx<nt_old; idx++)
		{
		  for (jdx = 0; jdx<ns_old; jdx++)
		    {
		      temp[nt_old*jdx+ idx] = projdbt[index+ns_old*idx+ jdx]; // transpose
		    }
		}
  
	      for (i = 0; i < ns*nt; i++) { proji[i]= 0.f; }
	      printf("view # %i\n", kdx);
	      // interpolate into CT geometry
	      ct_projections(temp, proji, betas_rad, kdx, ns, nt, ns_old, nt_old, offset_s, offset_t, offset_sold, offset_told, dt, ds, dod, dso, dsd, orbit); 

	      indexi = kdx*nt*ns;
	      for (i=0; i<nt*ns; i++)
		{
		  proj[i + indexi] = proji [i];
		}

	    }
	  free (proji); free(projdbt); free(temp);
	} // runlvl <= PREWEIGHT if ends



//////////////////////////////////TO DEBUG WITH DIFFERENT RUN LEVELS ///////////////////////////////////////

      if (runlvl == PREFILTER) 
	{
	  fread(proj, sizeof(double), ns*nt*na, fid);
	  printf("DONE READING INPUT!\n"); // Finish reading input
	}

      if (runlvl == BACKPROJ) 
	{
	  fread(proj, sizeof(double), ns*nt*na, fid);
	  printf("DONE READING INPUT!\n"); // Finish reading input
	}
      

//////////////////////////STARTING STANDARD FBP RECONSTRUCTION //////////////////////////////////////////
   // Weight input projections
      if (runlvl <= PREWEIGHT) 
	{
	  printf("Weighting projections ... ");
	  tac = clock();
	  fdk_weight( proj, ns, nt, na, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl);   
	  toc = clock();
	  tinc = (double)(toc-tac) / (double)CLOCKS_PER_SEC;
	  tspan += tinc;
	  printf("%lg sec\n", tinc);
	}



   // Filter input projections
      if (runlvl <= PREFILTER) 
	{
	  printf("Filtering projections ... ");
	  tac = clock();
	  fdk_filter(proj, window, dsd, dfs, ds, ns, nt, na);
	  toc = clock();
	  tinc = (double)(toc-tac) / (double)CLOCKS_PER_SEC;
	  tspan += tinc;
	  printf("%lg sec\n", tinc);   
	}

   

   // Perform back projection
      if (runlvl <= BACKPROJ) 
	{
	  printf("Reconstructing image .... ");
	  tac = clock();
	  cbct_back(  proj, ns, nt, na, betas_rad, offset_xyz, nx, nz, ny, dx, dz, -dy,
		      offset_s, offset_t, img_temp, ia_skip, dso, dfs, dsd, ds, dt, orbit);
	  toc = clock();
	  tinc = (double)(toc-tac) / (double)CLOCKS_PER_SEC;
	  tspan += tinc;
	  printf("%lg sec\n", tinc);

	  index = 0;
	  for (kdx = 0; kdx< nz; kdx++)
	    {
	      index = kdx*ny;   
	      for (idx=0; idx<nx; idx++)
		{
		  for (jdx = 0; jdx<ny; jdx++)
		    {
		      img_tr[index + idx*ny*nz + jdx] = img_temp[jdx*nx*nz + idx*nz+ kdx];  // 3D Transpose "permute.m"
		    }
		}
	    }

	  for (kdx = 0; kdx< nz; kdx++)
	    {
	      index = kdx*ny*nx;   
	      for (idx=0; idx<nx; idx++)
		{
		  for (jdx = 0; jdx<ny; jdx++)
		    {
		      img[index+nx*jdx + idx] = img_tr[index+ny*idx+ jdx]; // transpose
		    }
		}
	    }

	  for (idx = 0; idx < nx*ny*nz; idx++) { img[idx] = img[idx]*((double)180/orbit); }
	  
	} // runlvl <= BACKPROJ if ends

      // Write output image to file
      printf("Writing output image to file ... ");
      fid = fopen(outFile, "wb");
      if (fid == NULL) { printf("\nERROR: Could not open output file %s\n", outFile); }
      fwrite(img, sizeof(double), nx*ny*nz, fid);
      fclose(fid);
      printf("DONE!\n"); // Finish writing output

    
  
      toc = clock();                                                
      printf("Total processing time elapsed: %lg sec\n", tspan);   // total time include back-projection, weighting and filtering 
      tspan = (double)(toc-tic) / (double)CLOCKS_PER_SEC;          
      printf("Total execution time elapsed : %lg sec\n", tspan);   // total time for reading, interpolating, standard FBP and saving data
      free(proj);
      free(betas_rad);
      free(img);
      free(img_temp);
      free(img_tr);
   
    } // for (phntm = 0; phntm < f; phntm++) ends

  //fclose (fp);
  return 0;
}  // main ends 

