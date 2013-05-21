//============================================================
//  File:       jgp2pgm.cpp
//  Author:     John Gauch
// Modified Author: Lihang Li
//  Date:       May 2013
//============================================================


#include "jpg2pgm.h"

int main(int argc, char** argv)
{
    if(argc !=3)
    {
      printf("Usage: jpg2pgm infile outfile\n");
    
      return -1;
    }

    // Read input image
    im_short Image1;
    Image1.ReadJpg(argv[1]);

    // Write output image
    Image1.WriteBinary(argv[2]);
    return 0;
}

