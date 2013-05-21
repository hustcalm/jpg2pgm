//============================================================
//  File:       im_pipeline.h
//  Author:     John Gauch
//  Date:       Spring 2009
//============================================================

#ifndef _PIPELINE_
#define _PIPELINE_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "im_color.h"
#include "im_queue.h"

#define NUM_THREAD 10
#define NUM_QUEUE 10

class im_pipeline
{
 public:
   // Constructor functions
   im_pipeline();
   ~im_pipeline();

   // Input output functions
   void ReadJpg(int out_q, char *dir);
   void WriteJpg(int in_q, char *dir);

   // Imaging operations
   void Noise(int in_q, int out_q, float num);
   void Gaussian(int in_q, int out_q, float num);

 private:
     pthread_t Thread[NUM_THREAD];
   int ThreadCount;
   im_queue Queue[NUM_QUEUE];
   int QueueCount;
};

#endif
