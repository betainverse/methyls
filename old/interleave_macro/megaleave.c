/* megaleave.c --> from grad_sort_nd.c 
    
    Written by Lewis E. Kay on July 5 1992 to carry out time-domain 
            manipulations on Rance-type HSQC data to get pure-absn
            lineshape

    Modified by L.E.Kay on December 13 1992 to carry out time-domain
            manipulations on the ni2 dimension of a phase2,phase
            3D data set to get pure absn lineshape.

    Modified by L.E.Kay on April 11 1993 to carry out time-domain
             manipulations on 4D data sets acquired as d3,d2,d4,
             phase3,phase2,phase data sets.

    Modified by L.E.Kay on May 4 1993 to carry out time-domain
             manipulations on 2D data sets
	     
    Modified by M.R.Gryk on August 16 1993 to split and gradsort
     	     interleaved 2D data sets	   
	     
    Modified by M.R.Gryk on November 2 1998 to split more than 2 data
             sets and to flag for phase array position  
	    
    Modified by M.R.Gryk on May 15, 2001 to remove automatic gradsort
             operation for use with NMRPipe.

*/

#include <math.h>
#include <stdio.h>

#define MAXDATA 4097
#define MAXSTRING 20
#define MAXOUT 50             /* maximum number of interleaved spectra  */

  main(argc,argv)

     int argc;
     char *argv[];

   {
     
      FILE *fi, *fo[MAXOUT];
      char head[MAXOUT*2][50];
      char outname[MAXOUT][MAXSTRING];
      int data[MAXOUT*2][MAXDATA],
          sum[MAXOUT][MAXDATA],
          diff[MAXOUT][MAXDATA],
          phase[MAXOUT][MAXDATA],
          numblocks,jj,i,mg,ni,np,nsets, flag;


      if(argc !=7) {
        printf("\n TO BE USED WITH VARIAN DATA\n");
        printf("Input data set\n");
        printf("Output data prefix\n");
	printf("Number of data sets\n");
        printf("What is ni (integer)\n");
        printf("Total number of acq. points: np (REAL + IMAG)\n");
	printf("Array flag (0 or 1)\n"); 
	printf("  flag 0: phase is arrayed second, ie. phase, other\n");
	printf("  flag 1: phase is arrayed first,  ie. other, phase\n");
        exit(1);
      }

      fi = fopen(argv[1],"r");
      nsets = atoi(argv[3]);
      for (mg=0; mg<nsets; mg++) {
        sprintf(outname[mg], "%s_%d", argv[2], mg); 
        fo[mg] = fopen(outname[mg],"w");
      }
      ni = atoi(argv[4]);
      np = atoi(argv[5]);
      flag = atoi(argv[6]);


      /* read header */
      fread(&numblocks,1,4,fi);     /* first 4 bytes tells # blocks        */
      printf("%d blocks total, %d per expt\n", numblocks, numblocks/nsets);
      numblocks/=nsets;             /* divide this when split spectra      */
      fread(head[0],1,28,fi);       /* 28 more bytes in header             */
      for (mg=0; mg<nsets; mg++) {
        fwrite(&numblocks,1,4,fo[mg]);
        fwrite(head[0],1,28,fo[mg]);
      }

      for(jj=1; jj<=ni; jj++) {       /* read 1st 4 FID's: 1,3 are same exp  */
                                      /*  2,4 are same exp                   */
        for(mg=0; mg<(nsets*2); mg++) {
          fread(head[mg],1,28,fi);
          fread(&data[mg][1],4,np,fi);
	}
        
	
	for(mg=0; mg<nsets; mg++) {
	  if (flag == 0) {
             fwrite(head[mg],1,28,fo[mg]);
             fwrite(&data[mg][1],4,np,fo[mg]);
             fwrite(head[mg+nsets],1,28,fo[mg]);
             fwrite(&data[mg+nsets][1],4,np,fo[mg]);
	  } else {
	     fwrite(head[mg*2],1,28,fo[mg]);
             fwrite(&data[mg*2][1],4,np,fo[mg]);
             fwrite(head[mg*2+1],1,28,fo[mg]);
             fwrite(&data[mg*2+1][1],4,np,fo[mg]);
	  }
	}
      }


      fclose(fi); 
      for(mg=0; mg<nsets; mg++) fclose(fo[mg]);
      
  }
