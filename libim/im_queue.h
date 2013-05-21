//============================================================
//  File:       im_queue.h
//  Author:     John Gauch
//  Date:       Spring 2009
//============================================================

#ifndef _IM_QUEUE_
#define _IM_QUEUE_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "im_color.h"

#define MAX_LENGTH 10

struct im_node
{
   im_color Image;
   int Number;
   im_node *Next;
   im_node *Prev;
};

class im_queue
{
 public:
   // Constructor functions
   im_queue();
   ~im_queue();

   // Queue functions
   bool IsFull();
   bool IsEmpty();
   void WaitNotFull();
   void WaitNotEmpty();
   void Insert(int number, im_color & im);
   void Remove(int &number, im_color & im);
   int GetLength();

 private:
     pthread_mutex_t * Mutex;
   pthread_cond_t *NotFull;
   pthread_cond_t *NotEmpty;
   im_node *Tail;
   im_node *Head;
   int Length;
};

#endif
