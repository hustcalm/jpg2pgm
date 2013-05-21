//============================================================
//  File:       im_pipeline.cpp
//  Author:     John Gauch
//  Date:       Spring 2009
//============================================================

#include "im_pipeline.h"

//============================================================
im_pipeline::im_pipeline()
{
   ThreadCount = 0;
   QueueCount = 0;
}

//============================================================ 
im_pipeline::~im_pipeline()
{
   ThreadCount = 0;
   QueueCount = 0;
}

//============================================================
void *ReadJpg_Thread(void *data)
{
   // Unpack parameters for pthread
   struct Param
   {
      im_queue *out_q;
      char *dir;
   };
   Param *param = (Param *) data;

   // Get list of files
   char name[STRLEN];
   sprintf(name, "%s/files.txt", param->dir);
   char command[STRLEN];
   sprintf(command, "find %s -iname \"*jpg\" | sort > %s", param->dir, name);
   if (system(command) == -1) return NULL;

   // Process data
   int n = 0;
   FILE *fd = fopen(name, "r");
   while ((fd != NULL) && (fscanf(fd, "%s", name) != EOF))
   {
      im_color image;
      image.ReadJpg(name);
      param->out_q->WaitNotFull();
      param->out_q->Insert(n++, image);
      printf("ReadJpg %s\n", name);
   }
   return NULL;
}

//============================================================
void im_pipeline::ReadJpg(int out_q, char *dir)
{
   // Pack parameters for pthread
   struct Param
   {
      im_queue *out_q;
      char *dir;
   };
   Param *param = new Param;
   param->out_q = &(Queue[out_q]);
   param->dir = dir;

   // Create thread
   if (ThreadCount < NUM_THREAD)
      pthread_create(&Thread[ThreadCount++], NULL,
		     ReadJpg_Thread, (void *) param);
}

//============================================================
void *WriteJpg_Thread(void *data)
{
   // Unpack parameters for pthread
   struct Param
   {
      im_queue *in_q;
      char *dir;
   };
   Param *param = (Param *) data;

   // Create directory
   char command[STRLEN];
   sprintf(command, "mkdir -p %s", param->dir);
   if (system(command) == -1) return NULL;

   // Process data
   while (true)
   {
      int n;
      im_color image;
      param->in_q->WaitNotEmpty();
      param->in_q->Remove(n, image);
      char name[STRLEN];
      sprintf(name, "%s/%03d.jpg", param->dir, n);
      image.WriteJpg(name);
      printf("WriteJpg %d\n", n);
   }
   return NULL;
}

//============================================================
void im_pipeline::WriteJpg(int in_q, char *dir)
{
   // Pack parameters for pthread
   struct Param
   {
      im_queue *in_q;
      char *dir;
   };
   Param *param = new Param;
   param->in_q = &(Queue[in_q]);
   param->dir = dir;

   // Create thread
   if (ThreadCount < NUM_THREAD)
      pthread_create(&Thread[ThreadCount++], NULL,
		     WriteJpg_Thread, (void *) param);
}

//============================================================
void *Noise_Thread(void *data)
{
   // Unpack parameters for pthread
   struct Param
   {
      im_queue *in_q;
      im_queue *out_q;
      float num;
   };
   Param *param = (Param *) data;

   // Process data
   while (true)
   {
      int n;
      im_color image;
      param->in_q->WaitNotEmpty();
      param->in_q->Remove(n, image);
      image.NoiseImpulse('b', param->num);
      param->out_q->WaitNotFull();
      param->out_q->Insert(n, image);
      printf("Noise %d\n", n);
   }
   return NULL;
}

//============================================================
void im_pipeline::Noise(int in_q, int out_q, float num)
{
   // Pack parameters for pthread
   struct Param
   {
      im_queue *in_q;
      im_queue *out_q;
      float num;
   };
   Param *param = new Param;
   param->in_q = &(Queue[in_q]);
   param->out_q = &(Queue[out_q]);
   param->num = num;

   // Create thread
   if (ThreadCount < NUM_THREAD)
      pthread_create(&Thread[ThreadCount++], NULL,
		     Noise_Thread, (void *) param);
}

//============================================================
void *Gaussian_Thread(void *data)
{
   // Unpack parameters for pthread
   struct Param
   {
      im_queue *in_q;
      im_queue *out_q;
      float num;
   };
   Param *param = (Param *) data;

   // Process data
   while (true)
   {
      int n;
      im_color image;
      param->in_q->WaitNotEmpty();
      param->in_q->Remove(n, image);
      image.Gaussian(param->num);
      param->out_q->WaitNotFull();
      param->out_q->Insert(n, image);
      printf("Gaussian %d\n", n);
   }
   return NULL;
}

//============================================================
void im_pipeline::Gaussian(int in_q, int out_q, float num)
{
   // Pack parameters for pthread
   struct Param
   {
      im_queue *in_q;
      im_queue *out_q;
      float num;
   };
   Param *param = new Param;
   param->in_q = &(Queue[in_q]);
   param->out_q = &(Queue[out_q]);
   param->num = num;

   // Create thread
   if (ThreadCount < NUM_THREAD)
      pthread_create(&Thread[ThreadCount++], NULL,
		     Gaussian_Thread, (void *) param);
}
