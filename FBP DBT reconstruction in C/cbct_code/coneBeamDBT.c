//#*********************************************************************************************************************************
//#*********************************************************************************************************************************
//#
//#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
//#								
//# BASED ON THE C RECONSTRUCTION CODE DEVELOPED BY LEESER, MUKHERJEE AND BROCK (BMC RESEARCH NOTES 7.1, p. 582, 2014).
//# WHICH IN TURN IMPLEMENTS THE FDK RECONSTRUCTION ALGORITHM DEVELOPED BY FESSLER (http://web.eecs.umich.edu/~fessler/).
//#
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
//#							FDA DISCLAIMER
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
 This program takes a number of 2D projections and compute a 3D fron those.
 It has three steps : Weighting, Filtering and Back projection.
 Usage: running the makefile should work
 Input path is hardcoded in the source code.
 
 Copyright (C) 2012 James Brock
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "coneBeam.h"

//
// Cone beam CT back projection algorithm
//
void cbct_back(double * proj, int ns, int nt, int na, double * betas_rad, double * offset_xyz, int nx, int ny, int nz, double dx, double dy, double dz, double offset_s, double offset_t, double * img, int ia_skip, double dso, double dfs, double dsd, double ds, double dt) {
   int idx, jdx, kdx = 0;
   double offset_source = 0;
   double * source_zs = (double*)malloc(sizeof(double)*na);
   memset(source_zs, 0, sizeof(double)*na);
   int scale_dang = 1;
   int orbit = 360;
   int orbit_start = 0;
   clock_t tic, toc;
   double tspan = 0.0;

   
	
   

   double wx = (double(nx-1)/2.d) + offset_xyz[0];
   double wy = (double(ny-1)/2.d) + offset_xyz[1];
   double wz = (double(nz-1)/2.d) + offset_xyz[2];

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

 

   double * zc = (double*)malloc(sizeof(double)*nz);
   for (idx = 0; idx < nz; idx++) {
      zc[idx] = (double(idx) - wz) * dz;
   }
 
   double ws = (double(ns+1)/2.d) + offset_s;
   double wt = (double(nt+1)/2.d) + offset_t;

   // sdim = [ns+3 nt+3]; % trick: extra zeros saves indexing in loop
   int sdim[2] = {0, 0};
   sdim[0] = ns+3; sdim[1] = nt+3;
   double * proj1 = (double*)malloc(sizeof(double)*sdim[0]*sdim[1]);
   memset(proj1, 0, sizeof(double)*sdim[0]*sdim[1]);
   
   // for iz=1:nz
   for (idx = 0; idx < nz; idx++) {
   	int ia_min = 0;
      int ia_max = na-1;
      //	   printf("\tImage slice %d ... ", idx);
	   tic = clock();
	   // % loop over each projection angle
	   double * img2 = (double*)malloc(sizeof(double)*nx*ny);
	   memset(img2, 0, sizeof(double)*nx*ny);
		for (jdx = ia_min; jdx < ia_max; jdx += ia_skip) {
         double beta = betas_rad[jdx];

		   double * x_betas = (double*)malloc(sizeof(double)*nx*ny);
		   for (kdx = 0; kdx < nx*ny; kdx++) {
		      x_betas[kdx] = (xc[kdx]*cosf(beta)) + (yc[kdx]*sinf(beta));
		   }
		   double * y_betas = (double*)malloc(sizeof(double)*nx*ny);
         for (kdx = 0; kdx < nx*ny; kdx++) {
            y_betas[kdx] = dso - ((-1.d*xc[kdx]*sinf(beta)) + (yc[kdx]*cosf(beta)));
         }

		   // % detector indices
		   double * mag = (double*)malloc(sizeof(double)*nx*ny);
			if ((dsd==INFINITY) || (dso==INFINITY)) {
            for (kdx = 0; kdx < nx*ny; kdx++) {
               mag[kdx] = 1.d;
            }
		   } else {
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
			} else if (dfs == 0.d) {
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
		/*bs(bs < 1 | bs > ns) = ns+1; % trick for zero extrapolation in s
		bt = max(bt, 1); % implicit extrapolation in t (detector row)
		bt = min(bt, nt);
		*/
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

/*		  is(is == ns+1) = ns; % trick for zero extrapolation in s
		it(it == nt) = nt - 1; % trick for last row extrapolation in t
*/
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
		      wl[kdx] = 1.d - wr[kdx];
		   }
		   // upper weight
		   double * wu = (double*)malloc(sizeof(double)*nx*ny);
		   for (kdx = 0; kdx < nx*ny; kdx++) {
		      wu[kdx] = bt[kdx] - (double)it[kdx];
		   }
		   // lower weight
		   double * wd = (double*)malloc(sizeof(double)*nx*ny);
		   for (kdx = 0; kdx < nx*ny; kdx++) {
		      wd[kdx] = 1.d - wu[kdx];
		   }

		   unsigned int * ibad = (unsigned int*)malloc(sizeof(unsigned int)*nx*ny);
		   for (kdx = 0; kdx < nx*ny; kdx++) {
		      ibad[kdx] = (is[kdx] < 0) | (is[kdx] > ns) | (it[kdx] < 0) | (it[kdx] > nt);
		   }
		   for (kdx = 0; kdx < nx*ny; kdx++) {
		      if (ibad[kdx]) { is[kdx] = ns+1; it[kdx] = nt+1; } 
		   }

		   int mdx = 0;
		   int plane = ns*nt*jdx;
         for (kdx = 0; kdx < nt; kdx++) {
            int col = kdx*ns;
            int colx = (kdx+1)*sdim[0];
            for (mdx = 0; mdx < ns; mdx++) {
                proj1[colx+mdx+1] = proj[plane+col+mdx];
            }
         }
      
         
		   // vertical interpolation
			double * p1 = (double*)malloc(sizeof(double)*nx*ny);
			double * p2 = (double*)malloc(sizeof(double)*nx*ny);
			double * p0 = (double*)malloc(sizeof(double)*nx*ny);
			for (kdx = 0; kdx < nx*ny; kdx++) {
			   int is1 = is[kdx]; int is2 = is1 + 1;
			   int it1 = it[kdx]; int it2 = it1 + 1;
			   p1[kdx] = (wl[kdx] * proj1[(it1*sdim[0])+is1]) + (wr[kdx] * proj1[(it1*sdim[0])+is2]);
			   p2[kdx] = (wl[kdx] * proj1[(it2*sdim[0])+is1]) + (wr[kdx] * proj1[(it2*sdim[0])+is2]);
			   p0[kdx] = (wd[kdx] * p1[kdx]) + (wu[kdx] * p2[kdx]);
			}

			if ((dfs == INFINITY) || (dsd == INFINITY) || (dso == INFINITY)) {
			   for (kdx = 0; kdx < nx*ny; kdx++) {
			      p0[kdx] = p0[kdx] * mag[kdx] * mag[kdx];
			   }
			} else if (dfs == 0.d) {
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
        img[plane+kdx] = img2[kdx];
      }
      free(img2);
      
	   if (scale_dang) {
         for (kdx = 0; kdx < nx*ny; kdx++) { img[plane+kdx] =  0.1*img[plane+kdx] * ( (0.5f*(abs(orbit)*M_PI/180)) / (double(na)/double(ia_skip)) ); }
      }

         
      toc = clock();
      tspan = (double)(toc-tic) / (double)CLOCKS_PER_SEC;
      //printf("%f\n", tspan);
   } // end iz (slices)

   free(source_zs);
   free(xc); free(yc); free(zc);
   free(proj1);
   return;
}

//
// Feldkamp Weight algorithm
//
void fdk_weight(  double * proj, int ns, int nt, int na, double ds, double dt, double offset_s,
                  double offset_t, double dsd, double dso, double dfs, int w1cyl) {
   int idx, jdx, kdx, offset, index = 0;
   double tmp1, tmp2 = 0.d;

   // [ns nt na] = size(proj);
   // ss = ([-(ns-1)/2:(ns-1)/2]' - offset_s) * ds;
   // tt = ([-(nt-1)/2:(nt-1)/2]' - offset_t) * dt;
   // [ss tt] = ndgrid(ss, tt);
   double * ss = (double*)malloc(sizeof(double)*ns*nt);
   for (idx = 0; idx < ns; idx++) {
      tmp1 = ((idx - ((double(ns)-1.d) / 2.d)) - offset_s) * ds;
      for (jdx = 0; jdx < nt; jdx++) { ss[(jdx*ns)+idx] = tmp1; }
   }
   double * tt = (double*)malloc(sizeof(double)*ns*nt);
   for (idx = 0; idx < nt; idx++) {
      tmp2 = ((idx - ((double(nt)-1.d) / 2.d)) - offset_t) * dt;
      offset = idx*ns;
      for (jdx = 0; jdx < ns; jdx++) { tt[offset+jdx] = tmp2; }
   }
   
   double * ww = (double*)malloc(sizeof(double)*ns*nt);
   // if isinf(dfs) % flat
   if (dfs == INFINITY) {
      //ww1 = dso ./ sqrt(dsd^2 + ss.^2 + tt.^2);   
      for (idx = 0; idx < nt; idx++) {
         offset = idx*ns;
         for (jdx = 0; jdx < ns; jdx++) {
            index = offset+jdx;
            tmp1 = ss[index]*ss[index];
            tmp2 = tt[index]*tt[index];
            ww[offset+jdx] = dso / sqrt(dsd*dsd + tmp1 + tmp2);
         }
      }
	   // if ~w1cyl //% Jinyi Qi suggested not to do this for cylinder objects
		   // ww1 = ww1 .* sqrt(1 + (tt/dsd).^2); //% original version in book
	   // end
      if (!w1cyl) {
         for (idx = 0; idx < nt; idx++) {
            offset = idx*ns;
            for (jdx = 0; jdx < ns; jdx++) {
               index = offset + jdx;
               tmp1 = (tt[index] / dsd);
               ww[index] *= sqrt(1.d + (tmp1*tmp1));
            }
         }
      }
   // elseif dfs == 0 % arc
	} else if (dfs == 0.d) {
           if (!w1cyl){
	   // ww1 = (dso/dsd) * cos(ss ./ (dsd * sqrt(1 + (tt/dsd).^2)));
	   for (idx = 0; idx < nt; idx++) {
	      offset = idx*ns;
	      for (jdx = 0; jdx < ns; jdx++) {
	         index = offset + jdx;
	         tmp1 = (tt[index] / dsd);
	         ww[index] = (dso/dsd) * cos(ss[index] / (dsd * sqrt(1.d + (tmp1*tmp1))));
	      }
	   }
	}
	   // if w1cyl
		   // ww1 = ww1 ./ sqrt(1 + (tt/dsd).^2); % todo: new version
	   // end
	   else if (w1cyl) {
	      for (idx = 0; idx < nt; idx++) {
	         offset = idx*ns;
	         for (jdx = 0; jdx < ns; jdx++) {
	            index = offset + jdx;
	            tmp1 = (tt[index] / dsd);
	         //   ww[index] /= sqrt(1.d + (tmp1*tmp1));                                          // commented May 25
                    ww[index] = (dso/dsd) * cos(ss[index]/dsd)/sqrt(1.d + (tmp1*tmp1));           // edited May 25
                  
	         }
	      }
	   }
   // else
	} else {
	   printf("ERROR: Other configurations not implemented\n");
   // end
   }

   // for ia=1:na % same weighting for each view angle
	   // proj(:,:,ia) = proj(:,:,ia) .* ww1;
   // end
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
// if n/2 * ds / dsd > 0.9 * pi/2
//	   printm('angle is %g degrees: too large', rad2deg(n/2 * ds / dsd))
//	   error 'physically impossible arc geometry'
// end
   int idx = 0;
   double tmp = 0;
   // nn = [-(n/2):(n/2-1)]';
   int * nn = (int*)malloc(sizeof(int)*n);
   for (idx = 0; idx < n; idx++) {
      nn[idx] = idx - (n/2);
      // h = zeros(size(nn));
      // h(nn==0) = 1 / (4 * ds^2);
      if (nn[idx] == 0) {
         hval[idx] = 1.d / (4.d*ds*ds);
      } else {
         hval[idx] = 0.d;
      }
   }
   // odd = mod(nn,2) == 1;
   // h(odd) = -1 ./ (pi * dsd * sin(nn(odd) * ds / dsd)).^2;
   for (idx = 0; idx < n; idx++) {
      if ((nn[idx]%2) != 0) {
         tmp = sin(((double)(nn[idx]) * ds) / dsd) * dsd * M_PI;
         hval[idx] = -1.d / (tmp*tmp);
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
   // nn = [-(n/2):(n/2-1)]';
   int * nn = (int*)malloc(sizeof(int)*n);
   for (idx = 0; idx < n; idx++) {
      nn[idx] = idx - (n/2);
      // h = zeros(size(nn));
      hval[idx] = 0.d;
   }
   // h(n/2+1) = 1 / 4;
   hval[n/2] = 1.d / 4.d;
   // odd = mod(nn,2) == 1;
   // h(odd) = -1 ./ (pi * nn(odd)).^2;
   // h = h / ds^2;
   for (idx = 0; idx < n; idx++) {
      if ((abs(nn[idx])%2) == 1) {
         hval[idx] = -1.d / ((M_PI * nn[idx]) * (M_PI * nn[idx]));
      }
      hval[idx] /= (ds*ds);
   }
   free(nn);
   return;
}

//
//  function [h, nn] = fbp_ramp(type, n, ds, dsd)
//
void fbp_ramp(double * hval, char * type, int n, double ds, double dsd) {
   // switch type
   // case 'arc'
	if (!strcmp(type,"arc")) {
	   ramp_arc(hval, n, ds, dsd);
   // case 'flat'
   } else if (!strcmp(type,"flt")) {
	   ramp_flat(hval, n, ds);
   // otherwise
	} else {
	   printf("ERROR: bad fan type\n");
   // end
   }
   return;
}
   
//
// Feldkamp-Davis-Kress fan filter function
// 
void fdk_fan_filter(double * hval, char * type, int n, double ds, double dsd, char * wintype ) {
   const char * mFile = "./data/hval.txt";
   FILE * fid = NULL;
   int idx = 0;
   int offset = n/2;
   float * h_shift = (float*)malloc(sizeof(float)*n);
   fftwf_complex * fftc_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*((n/2)+1));
   fftwf_plan plan = fftwf_plan_dft_r2c_1d(n, h_shift, fftc_out, 0);

   fbp_ramp(hval, type, n, ds, dsd);

   // H = reale(fft(fftshift(h)));
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
   // if streq(window, 'ramp')
	if (!strcmp(wintype, "ramp")) {
	   // window = ones(n,1);
	   for (idx = 0; idx < n; idx++) { window[idx] = 1.d; }
   // elseif streq(window, 'hann')
   } else if (!strcmp(wintype,"hann")) {
	   //window = hann(n, 'periodic');
      for (idx = 0; idx < n; idx++) { window[idx] = 0.5f * (1.d + cos((2.d*M_PI*idx)/(n-1))); }
   // else
	} else {
	   // window = fftshift(fbp2_window(n, window));
      printf("ERROR: Unknown window type!\n");
   // end
   }
   
  

   // H = H .* fftshift(window);
   for (idx = 0; idx < n; idx++) {
      if (idx < (n/2)) {
         h_shift[idx] = window[idx+offset];
      } else {
         h_shift[idx] = window[idx-offset];
      }
      hval[idx] *= h_shift[idx];
   }
/*
   printf("Writing filter to file ... ");
   fid = fopen(mFile, "wb");
   if (fid == NULL) { printf("\nERROR: Could not open output file %s\n", mFile); }
   fwrite(hval, sizeof(double), n, fid);
   fclose(fid);
*/
   fftwf_destroy_plan(plan);
   fftwf_free(fftc_out);
   fftwf_free(h_shift);
   return;
}

//
// Feldkamp-Davis-Kress filter function
//
void fdk_filter(  double * proj, char * window, double dsd, double dfs, double ds,
                  int ns, int nt, int na) {
   int idx, jdx, kdx, offset, index = 0;
   int npad = 1<<(int)ceil(log2((2*ns)-1));
   //printf ("npad %i \n", npad);
   float * proj_pad = (float*)malloc(sizeof(float)*npad);
   for (idx = 0; idx < npad; idx++) { proj_pad[idx] = 0.d; }
   fftwf_complex * fftout = (fftwf_complex*)fftwf_malloc(sizeof(fftw_complex)*npad);
   fftwf_plan planF = fftwf_plan_dft_r2c_1d(npad, proj_pad, fftout, 0);
   fftwf_plan planR = fftwf_plan_dft_c2r_1d(npad, fftout, proj_pad, 0);

   double * hval  = (double*)malloc(sizeof(double)*npad);
   for (idx = 0; idx < npad; idx++) { hval[idx] = 0.d; }
   // if isinf(dsd) % parallel-beam
   if (dsd == INFINITY) {
	   // H = fdk_fan_filter('flat', npadh, ds, [], window); % [nb 1]
      fdk_fan_filter(hval, "flt", npad, ds, 0.d, window);
   // elseif isinf(dfs) % flat
	} else if (dfs == INFINITY) {
	   // H = fdk_fan_filter('flat', npadh, ds, [], window); % [nb 1]
      fdk_fan_filter(hval, "flt", npad, ds, 0.d, window);
   // elseif dfs == 0 % arc
	} else if (dfs == 0.d) {
	   // H = fdk_fan_filter('arc', npadh, ds, dsd, window);
      fdk_fan_filter(hval, "arc", npad, ds, dsd, window);
   // end
   }
   // H = ds * H; % differential for discrete-space convolution vs integral
   for (idx = 0; idx < npad; idx++) { hval[idx] *= ds; }

   // proj = ifft_sym( fft(proj, npadh, 1) .* repmat(H, [1 nt na]) );
   // proj = proj(1:ns,:,:); % back to original unpadded size
   for (idx = 0; idx < na; idx++) {
      offset = idx*ns*nt;
      for (jdx = 0; jdx < nt; jdx++) {
         index = offset+(jdx*ns);
         for (kdx = 0; kdx < ns; kdx++) { proj_pad[kdx] = proj[index+kdx]; }
         fftwf_execute(planF);
         /*for (kdx = 0; kdx < (npad/2); kdx++) {
            fftout[kdx+(npad/2)+1][0] = fftout[(npad/2)-1-kdx][0];
            fftout[kdx+(npad/2)+1][1] = fftout[(npad/2)-1-kdx][1] * -1.d;   
         }*/
         for (kdx = 0; kdx < npad; kdx++) {
            fftout[kdx][0] *= hval[kdx];
            fftout[kdx][1] *= hval[kdx];
         }
         fftwf_execute(planR);
         for (kdx = 0; kdx < npad; kdx++) { proj_pad[kdx] /=  (double)(npad); }
         for (kdx = 0; kdx < ns; kdx++) { proj[index+kdx] = proj_pad[kdx]; }
	/* for (kdx = 0; kdx < ns; kdx++) { 
         if (kdx<4){proj[index+kdx] = proj_pad[4]; }
         else if (kdx>ns-5){proj[index+kdx] = proj_pad[ns-5];}
	}*/
      }
   }

   fftwf_destroy_plan(planF);
   fftwf_destroy_plan(planR);
   fftwf_free(fftout);
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
void dbt_to_ct(  double * projv, double * proji, double * betas_rad, int i, int ns, int nt, int ns_old, int nt_old, double offset_s, double offset_t, double offset_sold, double offset_told, double dt, double ds, double dod, double dso){
	
	double beta, F1, F2, F3, F4, ang, acos, asin, ws, wt, ws_old, wt_old, xs, ys, zs = 0;
        int idx, jdx, nrows, ncols, flagds, flagdt, dis, djs, dit, djt;
  	FILE * fid = NULL;
        beta = betas_rad[i];
        ang = (beta*3.141592653589793)/180;
        //ang  = (-15*3.141592653589793)/180;
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
// Main CBCT
//
int main(int argc, char *argv[]) {
   // CBCT Variables
   int index, indexi, ns_old, nt_old, ns, nt, nviews, na, nx, ny, nz = 0;
   int i, j, x, z, a, b, c, jdx, kdx;
   double ds, dt, dx, dy, dz, dso, dfs, dod, dsd = 0.d;
   double offset_s, offset_t = 0.d;
   double offset_sold, offset_told = 0.d;
   int ia_skip, w1cyl, idx = 0;
   double offset_xyz[3] = {0.d, 0.d, 0.d};
   char window[5] = {'r', 'a', 'm', 'p', '\0'};
   double *proji, *proj, *projdbt, *temp, *betas_rad, *img, *img_cmp = NULL;
   double *angles;
   const char * oFile = "./data/projections_interpolatedf.dat";
   FILE * fid = NULL;
   clock_t tic, tac, toc = 0;
   double tinc, tspan = 0.0;
   char inFile[MAXSTRLEN], mFile[MAXSTRLEN];
   // uint phant = 4; // choose the phantom : use 8 commented June 13
   uint phant = 0;

   // Start timer
   tic = clock();

   // Read input variables from file based on run level
   printf("Reading input data ... \n");

   /* sprintf(inFile, "../data/mouse/mouse_code_params.dat"); */
   /* fid = fopen(inFile, "rb"); */
   /* if (fid == NULL) { printf("\nERROR: Could not open pre-weight file %s\n", inFile); } */

  if (runlvl == PREWEIGHT) {
      sprintf(inFile, "./data/phantom%d/cc/phantom_dbt%d_allint.dat", phant, phant);
      fid = fopen(inFile, "rb");
      if (fid == NULL) { printf("\nERROR: Could not open pre-weight file %s\n", inFile); }
   } else if (runlvl == PREFILTER) {
      sprintf(inFile, "./data/phantom%d/cc/phantom%d_wtd.dat", phant, phant);
      fid = fopen(inFile, "rb");
      if (fid == NULL) { printf("\nERROR: Could not open pre-filter file %s\n", inFile); }
   } else if (runlvl == BACKPROJ) {
      sprintf(inFile, "./data/phantom%d/cc/phantom%d_flt.dat", phant, phant);
      fid = fopen(inFile, "rb");
      if (fid == NULL) { printf("\nERROR: Could not open back projection file %s\n", inFile); }
   }

/////////////edited May 19 //////////////////////

 // na = 360;                       // the number of 2D projections
 // ns = 64;                        // Commented Jun 13. 
 //  nt = 60;   
    
    na = 7;
    nviews = 7;                       // the number of 2D projections
    ns_old = 614;                        
    nt_old = 169;
    ns = 696; 
    nt = 356;    

   nx = 64;
   ny = 60;
   nz = 50;
 //  ds = 1.600;                // Commented Jun 13.
 //  dt = 1.600;

   ds = 0.04;
   dt = -0.04;   

   dx = 0.7812500;
   dy = -0.7812500;
   dz = 0.7812500;
  // dso = 54.1;                // Commented Jun 13.
   dod  = 4.5;
   dso  = 60.5;
   dfs = 0;
   dsd = 94.900;
  // offset_s = 0.25;           // Commented Jun 13. 
   offset_s = 0.0;
   offset_t = 0.0;
   
   offset_sold = 0.0;
   offset_told = -85.0;
   
   ia_skip = 1;
   w1cyl   = 1;
   offset_xyz[1] = 0;
   offset_xyz[2] = 0;
   offset_xyz[3] = 0;

   projdbt = (double*)malloc(sizeof(double)*ns_old*nt_old*na);
   temp = (double*)malloc(sizeof(double)*ns_old*nt_old);
   proji = (double*)malloc(sizeof(double)*ns*nt);
   proj = (double*)malloc(sizeof(double)*ns*nt*na);
   fread(projdbt, sizeof(double), ns_old*nt_old*na, fid);
   betas_rad = (double*)malloc(sizeof(double)*na);
   fread(betas_rad, sizeof(double), na, fid);
   fclose(fid);  

   for (x=0; x<na; x++)
   printf ("angles %f x %i \n", betas_rad[x], x+1);
   printf("DONE READING INPUT!\n"); // Finish reading input
   // Initialize output image
   img = (double*)malloc(sizeof(double)*nx*ny*nz);
   for (idx = 0; idx < nx*ny*nz; idx++) { img[idx] = 0.d; }
   for (idx = 0; idx < na*ns*nt; idx++) { proj[idx]= 0.d; }

   for (kdx = 0; kdx< nviews; kdx++){
        index = kdx*nt_old*ns_old;   
       for (idx=0; idx<nt_old; idx++){
        for (jdx = 0; jdx<ns_old; jdx++){
         temp[nt_old*jdx+ idx] = projdbt[index+ns_old*idx+ jdx]; // transpose
        }
       }
  
       for (i = 0; i < ns*nt; i++) { proji[i]= 0.d; }
       
       // interpolate into CT geometry
       dbt_to_ct(temp, proji, betas_rad, kdx, ns, nt, ns_old, nt_old, offset_s, offset_t, offset_sold, offset_told, dt, ds, dod, dso); 
       indexi = kdx*nt*ns;
       for (i=0; i<nt*ns; i++){
          proj[i + indexi] = proji [i];
      }
    }
  free (proji); free(projdbt); free(temp);
   // Read in expected (Matlab) output
/*   sprintf(mFile, "../data/phantom%d/cc/phantom%d_rcn.dat", phant, phant);
   fid = fopen(mFile, "rb");
   if (fid == NULL) { printf("\nERROR: Could not open Matlab Image file %s\n", mFile); }
   fread(img_cmp, sizeof(double), nx*ny*nz, fid);
   fclose(fid);*/

/*   // Commented June 13 to Test DBT pre-processing. Uncomment till line :  printf("Total execution time elapsed : %lg sec\n", tspan);
 
   // Weight input projections
   if (runlvl <= PREWEIGHT) {
      printf("Weighting projections ... ");
      tac = clock();
      fdk_weight( proj, ns, nt, na, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl);   
      toc = clock();
      tinc = (double)(toc-tac) / (double)CLOCKS_PER_SEC;
      tspan += tinc;
      printf("%lg sec\n", tinc);
   }
   // Filter input projections
   if (runlvl <= PREFILTER) {
      printf("Filtering projections ... ");
      tac = clock();
      fdk_filter(proj, window, dsd, dfs, ds, ns, nt, na);
      toc = clock();
      tinc = (double)(toc-tac) / (double)CLOCKS_PER_SEC;
      tspan += tinc;
      printf("%lg sec\n", tinc);   
   }
   // Perform back projection
   if (runlvl <= BACKPROJ) {
      printf("Reconstructing image .... ");
      tac = clock();
      cbct_back(  proj, ns, nt, na, betas_rad, offset_xyz, nx, ny, nz, dx, dy, dz,
                  offset_s, offset_t, img, ia_skip, dso, dfs, dsd, ds, dt);
      toc = clock();
      tinc = (double)(toc-tac) / (double)CLOCKS_PER_SEC;
      tspan += tinc;
      printf("%lg sec\n", tinc);
   }

   // Write output image to file
   printf("Writing output image to file ... ");
   fid = fopen(oFile, "wb");
   if (fid == NULL) { printf("\nERROR: Could not open output file %s\n", oFile); }
   //fwrite(&nx, sizeof(int), 1, fid);
   //fwrite(&ny, sizeof(int), 1, fid);
   //fwrite(&nz, sizeof(int), 1, fid);
   fwrite(proj, sizeof(double), na*ns*nt, fid);
   fclose(fid);
   printf("DONE!\n"); // Finish writing output

//   toc = clock();                                                // Commented June 13
//   printf("Total processing time elapsed: %lg sec\n", tspan);    // Commented June 13
//   tspan = (double)(toc-tic) / (double)CLOCKS_PER_SEC;           // Commented June 13
//   printf("Total execution time elapsed : %lg sec\n", tspan);    // Commented June 13 
*/

 

   printf("Writing output image to file ... ");
   fid = fopen(oFile, "wb");
   if (fid == NULL) { printf("\nERROR: Could not open output file %s\n", oFile); }
   fwrite(proj, sizeof(double), ns*nt*na, fid);
   fclose(fid);
   printf("DONE!\n"); // Finish writing output


   free(proj);
   free(betas_rad);
   free(img);   

   //free(img_cmp);
   return 0;
}   

